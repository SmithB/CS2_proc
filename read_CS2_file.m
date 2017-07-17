function [D_POCA, D_swath, orb, out_table] = pick_and_run(L1_filename_full, hemisphere, bursts, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       PROCESS 20Hz SIN L1B DATA
%                           Using pick_slope.m
%
%              Find L1b retracking point and 3 ambiguities
%
% DEFINITIONS
%
% -----from L1 data:------
%
% CP.wave_watts    %waveforms corrected to bin x power
% CP.coherence     %coherence
% CP.phase         %phase difference (radians)
% CP.retrack_bin   %retracked waveform bin of POCA
% CP.retrack_PD    %retracked phase difference from POCA
% CP.range         %ranges (m)
% CP.inf_angle     %inferred angles from Z to POCA
%
% CP.xv_X
% CP.yv_X    %components of X (nadir vector)
% CP.zv_X
%
% CP.xv_V
% CP.yv_V    %components of V (velocity vector)
% CP.zv_V
%
% CP.xv_Y
% CP.yv_Y    %components of Y (ellipsoidal velocity vector)
% CP.zv_Y
%
% CP.xv_Z
% CP.yv_Z    %components of Z (range vector before angle correction)
% CP.zv_Z
%
% CP.xv_R_minus
% CP.yv_R_minus    %components of R (range vector, -1 ambiguity)
% CP.zv_R_minus
%
% CP.xv_R_0
% CP.yv_R_0        %components of R (range vector, 0 ambiguity)
% CP.zv_R_0
%
% CP.xv_R_plus
% CP.yv_R_plus     %components of R (range vector, +1 ambiguity)
% CP.zv_R_plus
%
% CP.xps_sc      %x coordinate of spacecraft (polar sterographic)
% CP.yps_sc      %y coordinate of spacecraft (polar sterographic)
% CP.h_sc        %height of  (m)
%
% CP.xps_minus   %x coordinate of new point (polar sterographic, -1 ambiguity)
% CP.yps_minus   %y coordinate of new point(polar sterographic, -1 ambiguity)
% CP.h_minus     %height of new point (m, -1 ambiguity)
%
% CP.xps_0       %x coordinate of new point (polar sterographic, 0 ambiguity)
% CP.yps_0       %y coordinate of new point(polar sterographic, 0 ambiguity)
% CP.h_0         %height of new point (m, 0 ambiguity)
%
% CP.xps_plus    %x coordinate of new point (polar sterographic, +1 ambiguity)
% CP.yps_plus    %y coordinate of new point(polar sterographic, +1 ambiguity)
% CP.h_plus      %height of new point (m, +1 ambiguity)
%
% CP.bad_file    %file not on land

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[RefEll.a,RefEll.b,RefEll.e2,RefEll.finv]=refell('WGS84');
out_table=[];

[HDR_L1b, L1b]=Cryo_L1b_read(L1_filename_full);


%----added by alex


orb.cycle = HDR_L1b.CYCLE; %Cycle number
orb.rel_orbit = HDR_L1b.REL_ORBIT; %Relative orbit number
orb.abs_orbit = HDR_L1b.ABS_ORBIT; %Absolute orbit number
CP.lat_sc = L1b.GEO.LAT(:);
CP.lon_sc = L1b.GEO.LON(:);

Nsamps=size(L1b.SIN.data,1);
%convert 20 Hz waveform from counts to power
scale_a = repmat(L1b.SIN.echo_scaling,[1 1 Nsamps]);
scale_a = permute(scale_a,[3 1 2]);
scale_b = repmat(L1b.SIN.echo_scale_power,[1 1 Nsamps]);
scale_b = permute(scale_b,[3 1 2]);
CP.wave_watts = L1b.SIN.data.*((scale_a*1e-9).*(2.^scale_b));
CP.coherence = L1b.SIN.coherence;
CP.phase = L1b.SIN.phase_difference./1e6;

[i, j, k] = size(CP.wave_watts);
CP.wave_watts = reshape(CP.wave_watts, i, j*k);
CP.coherence = reshape(CP.coherence, i, j*k);
CP.phase = reshape(CP.phase, i, j*k);
CP.noise_dB=L1b.MEA.noise_power(:); 
CP.noise_dB(CP.noise_dB==0)=NaN;
CP.sec_number=reshape(repmat((1:k), [j, 1]), [1, j*k]);  %index into 1-hz arrays for each burst

if ~exist('bursts','var')
    bursts=[];
end

%----------------------determine retracking point-------------------------%
% pick the slope points
D_POCA = pick_slope(CP.wave_watts, CP.coherence, CP.phase, bursts, CP.noise_dB, L1b);
 
% ------------------now run the location routine
D_POCA=locate_CS_returns(D_POCA, L1b, CP, hemisphere);
if ~isempty(D_POCA.h)
    if ~exist('params','var'); 
        D_POCA=[]; orb_table=[]; orb=[]; D_swath=[]; return; 
    end
    if ~isstruct(params.DEM);
        this_DEM=read_geotif_xy(params.DEM, range(real(D_POCA.xPS(:)))+[-90 90], range(imag(D_POCA.xPS(:)))+[-90 90]);
    else
        this_DEM=params.DEM;
    end
    D_POCA=select_swath_ambiguity(D_POCA, this_DEM, {'xPS','h','time','burst','ambiguity', 'phase','power','coherence', 'samp','ret_count','range_surf'});
end
if isempty(D_POCA); D_swath=[]; orb=[]; return; end

% 
D_swath=proc_CS2_swath(CP.wave_watts, CP.coherence, CP.phase, bursts);
D_swath=locate_CS_returns(D_swath, L1b, CP, hemisphere);
if ~isempty(D_swath.h)
    if hemisphere==1
        DEM=read_geotif_xy('/Volumes/insar2/gmap/gimp/gimp_90m_cubic.tif', range(real(D_swath.xPS(:)))+[-90 90], range(imag(D_swath.xPS(:)))+[-90 90]);
    else
        if isstruct(params.DEM); 
            DEM=params.DEM;
        elseif exist(params.DEM,'file')
            DEM=read_geotif_xy(params.DEM, range(real(D_swath.xPS(:)))+[-90 90], range(imag(D_swath.xPS(:)))+[-90 90]);
        end
    end
    
    D_swath.R_offnadir=abs(D_swath.xPS-repmat(D_swath.xps_sc+1i*D_swath.yps_sc, [1, 3]));    
    D_swath=select_swath_ambiguity(D_swath, DEM, {'xPS','h','time','burst','seg_ind','ambiguity', 'phase','power','coherence', 'samp','R_offnadir','range_surf'});
    
    if ~isempty(D_swath);
        % loop over the POCA points, find the closest POCA point (in y) to
        % each swath point, record the y difference
        % find the closest POCA point (in range), record the difference
        D_swath.R_POCA=NaN(size(D_swath.time));
        D_swath.dRange_POCA=NaN(size(D_swath.time));
        D_swath.POCA_flag=false(size(D_swath.time));
        uB=unique(D_POCA.burst);
        for kB=1:length(uB)
            iB=find(D_POCA.burst==uB(kB));
            iS=(D_swath.burst==uB(kB));
            this_R_POCA=NaN(sum(iS),1);
            this_dRange_POCA=NaN(sum(iS),1);
            for kkB=1:length(iB);
                this_R_POCA=min(this_R_POCA, abs(D_POCA.xPS(iB(kkB))-D_swath.xPS(iS)));
                this_dRange_POCA=min(this_dRange_POCA, abs(D_POCA.range_surf(iB(kkB))-D_swath.range_surf(iS)));
            end
            D_swath.R_POCA(iS)=this_R_POCA;
            D_swath.dRange_POCA(iS)=this_dRange_POCA;
            D_swath.POCA_flag(iS)=D_swath.range_surf(iS) < min(D_POCA.range_surf(iB));
        end
    end
else
    D_swath=[];
end
% write out flags, etc
flag_vals=false(size(L1b.SIN.data,2), size(L1b.SIN.data,3));
flag_loc=define_error_flags;

f=fieldnames(flag_loc.GEO.MCD_FLAG);
for kf=1:length(f);
    flag_vals=flag_vals | L1b.GEO.MCD_FLAG.(f{kf});
    %figure(kf); clf; imagesc(L1b.GEO.MCD_FLAG.(f{kf})); title(f{kf}); colorbar
end
f=fieldnames(flag_loc.GEO.INS_CFG);
for kf=1:length(f);
    flag_vals=flag_vals | L1b.GEO.INS_CFG.(f{kf});
end

f=fieldnames(flag_loc.SIN.FLAG);
for kf=1:length(f);
    flag_vals=flag_vals | L1b.SIN.FLAG.(f{kf});
end

 

D_POCA.error_composite=flag_vals(D_POCA.burst);
AD=sign(diff(CP.lat_sc)); AD(end+1)=AD(end); 
out_table.burst=1:length(flag_vals(:));
out_table.error_composite=flag_vals(:);
out_table.AD=AD(:);
D_POCA.AD=AD(D_POCA.burst); 

function flag_loc=define_error_flags

flag_loc=[];

flag_loc.GEO=struct('MCD_FLAG',[],'INS_CFG', []);
flag_loc.GEO.MCD_FLAG=struct(...
    'Block_Degraded', [], ...
    'Blank_Block', [], ...
    'Datation_Degraded', [], ...
    'Orbit_Propag_Err', [], ...
    'Orbit_Discontinuity', [], ...
    'Echo_Saturation', [], ...
    'Other_Echo_Err', [], ...
    'Rx1_Err_SARin', [], ...
    'Rx2_Err_SARin', [], ...
    'AGC_Incon', [], ...
    'CAL1_Corr_Miss', [], ...
    'CAL1_Corr_IPF', [], ...
    'DORIS_USO_Corr', [], ...
    'Complex_CAL1_Corr_IPF', [], ...
    'TRK_ECHO_Err', [], ...
    'RX1_ECHO_Err', [], ...
    'RX2_ECHO_Err', [], ...
    'CAL2_Corr_Miss', [], ...
    'CAL2_Corr_IPF', [], ...
    'Power_Scaling_Err', [], ...
    'Att_Corr_Miss', []);

flag_loc.GEO.INS_CFG=struct('Echo_Satur_Err', [], 'Cycle_Report_Err', []);

flag_loc.SIN=struct('FLAG', struct(...
    'Multilook_Incomplete', [], ...
    'Beam_Angle_Steering_Err', []));



