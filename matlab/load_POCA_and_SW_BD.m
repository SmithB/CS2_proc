function [D_POCA, D_swath, orb, out_table] = load_POCA_and_SW_BD(L1_filename_full, hemisphere, bursts, params)

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
D_POCA=[]; D_swath=[]; orb=[]; orb_table=[];
[RefEll.a,RefEll.b,RefEll.e2,RefEll.finv]=refell('WGS84');
out_table=[];

[D, orb]=read_CS2_BD(L1_filename_full);

CP.lat_sc = D.lat_20_ku;
CP.lon_sc = D.lon_20_ku;

% scaling for power is done in read_CS2_BaselineD
CP.wave_watts = D.power;
CP.coherence = D.coherence_waveform_20_ku;
CP.phase = D.ph_diff_waveform_20_ku;

CP.noise_dB=D.noise_power_20_ku(:); 
CP.noise_dB(CP.noise_dB==0)=NaN;
% sec_number may be needed.  It comes from ind_meas_1hz_20_ku

if ~exist('bursts','var')
    bursts=[];
end

%----------------------determine retracking point-------------------------%
% pick the slope points
D_POCA = pick_slope_BD(CP.wave_watts, CP.coherence, CP.phase, bursts, CP.noise_dB);
 
% ------------------now run the location routine
D_POCA=locate_CS_returns_BD(D_POCA, D, CP, hemisphere, params);
D_POCA.error_composite=D.flag_mcd_20_ku(D_POCA.burst);

if ~isempty(D_POCA.h)
    if ~exist('params','var') 
        D_POCA=[];  orb=[]; D_swath=[]; return; 
    end
    D_POCA=select_swath_ambiguity_BD(D_POCA,  {'xPS','h','time','burst','ambiguity', 'phase','power','coherence', 'samp','ret_count','range_surf', 'DEM', 'geoid'});
end
if isempty(D_POCA); D_swath=[]; orb=[]; return; end

% locate the CS points
D_swath=proc_CS2_swath_BD(CP.wave_watts, CP.coherence, CP.phase, bursts);
D_swath=locate_CS_returns_BD(D_swath, D, CP, hemisphere, params);

if ~isempty(D_swath.h)   
    D_swath.R_offnadir=abs(D_swath.xPS-repmat(D_swath.xps_sc+1i*D_swath.yps_sc, [1, 3]));    
    D_swath=select_swath_ambiguity_BD(D_swath, {'xPS','h','time','burst','seg_ind','ambiguity', 'phase','power','coherence', 'samp','R_offnadir','range_surf', 'DEM','geoid'});
    
    if ~isempty(D_swath) && ~isempty(D_swath.time)
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
            for kkB=1:length(iB)
                this_R_POCA=min(this_R_POCA, abs(D_POCA.xPS(iB(kkB))-D_swath.xPS(iS)));
                this_dRange_POCA=min(this_dRange_POCA, abs(D_POCA.range_surf(iB(kkB))-D_swath.range_surf(iS)));
            end
            D_swath.R_POCA(iS)=this_R_POCA;
            D_swath.dRange_POCA(iS)=this_dRange_POCA;
            D_swath.POCA_flag(iS)=D_swath.range_surf(iS) < min(D_POCA.range_surf(iB));
        end
    end
    D_swath.error_composite=D.flag_mcd_20_ku(D_swath.burst);

else
    D_swath=[];
end

