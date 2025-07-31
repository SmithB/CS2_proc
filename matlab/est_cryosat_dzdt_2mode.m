function varargout=est_cryosat_dzdt_2mode(varargin)

if nargout>0
    [varargout{1:nargout}]=feval( varargin{:});
else
    feval(varargin{:});
end

%----------------------------------------------------------------------------
function [XG, region_params]=define_region(region_name, delta, data_file)

INDEX=index_point_data_h5('read_index', data_file);
XG=unique(round_to(INDEX(1).bin_x+1i*INDEX(1).bin_y, delta));
region_params=[];
switch lower(region_name)
    case('amundsen')
        delta=2.5e4;
        W=delta+2e4;
        temp=read_gmt_vectors('/Volumes/ice1/ben/Cryosat/clip_files/amundsen_fit_poly_ll.gmt');
        amundsen_outline=ll2ps(temp.ll(:,2), temp.ll(:,1));
        %amundsen_outline=temp.ll(:, 1:2)*[1; i];
        XG=XG(in_polygon(amundsen_outline, XG));
    case('thw_lakes')
        XG=define_region('amundsen', delta, data_file);
        XR=[-1450000    -1350000]+[-2e5 2e5];
        YR= [-500000     -400000]+[-2e5 2e5];
        XG=XG(real(XG) >= XR(1) & real(XG) <= XR(2) & imag(XG) >= YR(1) & imag(XG)<= YR(2));
        region_params.E_RMS_d2z0dx2=2e-7;  % marginally strong enough to suppress small peaks
        region_params.E_RMS_d3zdx2dt=6e-8;
        region_params.E_RMS_d2zdt2=1;
        region_params.power_fit=false;
    case('thw_lakes_coarse');
        [XG, region_params]=define_region('thw_lakes',delta, data_file);
        region_params.E_RMS_d2z0dx2=region_params.E_RMS_d2z0dx2/4;
        region_params.E_RMS_d3zdx2dt=region_params.E_RMS_d3zdx2dt/4;
    case('pig')
        XG=define_region('amundsen');
        XR=[ -1720000    -1370000 ]+[-1 1]*2.5e4;
        YR= [ -330000      -60000 ]+[-1 1]*2.5e4;
        XG=XG(real(XG) >= XR(1) & real(XG) <= XR(2) & imag(XG) >= YR(1) & imag(XG)<= YR(2));
    case('belling')
        temp=read_gmt_vectors('/Volumes/ice1/ben/Cryosat/clip_files/belling_fit_poly_ll.gmt');
        belling_outline=ll2ps(temp.ll(:,2), temp.ll(:,1));
        %amundsen_outline=temp.ll(:, 1:2)*[1; i];
        XG=XG(in_polygon(belling_outline, XG));
        region_params.E_RMS_d2z0dx2=2e-7/10;
        region_params.E_RMS_d3zdx2dt=6e-8/10;
        region_params.E_RMS_d2zdt2=0.2;
    case('larsenc')
        %???
        temp=read_gmt_vectors('/Volumes/ice1/ben/Cryosat/clip_files/LarsenC_fit_poly_ll.gmt');
        this_outline=ll2ps(temp.ll(:,2), temp.ll(:,1));
        XG=XG(in_polygon(this_outline, XG));
        region_params.E_RMS_d2z0dx2=2e-7;
        region_params.E_RMS_d3zdx2dt=6e-8/10;
        region_params.E_RMS_d2zdt2=0.2;
    case('test')
        
        if 1% Jako
            XG=(-1.8000e+02 - 2.2700e+03i)*1e3;
            region_params.E_RMS_d2z0dx2=2e-6;
            region_params.E_RMS_d3zdx2dt=6e-8;
            region_params.E_RMS_d2zdt2=1;
        end
        
        if 0
            % a charming spot on the Abbott, where if editing goes wrong,
            % the curvature at the edge of Farwell Island leads to a
            % progressive deletion of data until the island expands to fill
            % the ice shelf.
            %XG=-1.4110e+06+1i*-4.3033e+05;
            XG=(-1825-1i*50)*1e3;
            region_params.E_RMS_d2z0dx2=2e-7/10;
            region_params.E_RMS_d3zdx2dt=6e-8/10;
            region_params.E_RMS_d2zdt2=0.2;
        end
 
         if 0
            % a charming spot on the Abbott, where if editing goes wrong,
            % the curvature at the edge of Farwell Island leads to a
            % progressive deletion of data until the island expands to fill
            % the ice shelf.  RELAXED PARAMETERS
            %XG=-1.4110e+06+1i*-4.3033e+05;
            XG=(-1825-1i*50)*1e3;
            region_params.E_RMS_d2z0dx2=2e-6;
            region_params.E_RMS_d3zdx2dt=6e-8;
            region_params.E_RMS_d2zdt2=0.5;
         end
         if 0
             temp=read_gmt_vectors('/Volumes/ice1/ben/Cryosat/clip_files/belling_fit_poly_ll.gmt');
             belling_outline=ll2ps(temp.ll(:,2), temp.ll(:,1));
             temp=read_gmt_vectors('/Volumes/ice1/ben/Cryosat/clip_files/Clip_amundsen.gmt');
             amundsen_outline=(temp.ll(:,1)+1i*temp.ll(:,2));
             %amundsen_outline=temp.ll(:, 1:2)*[1; i];
             XG=XG(in_polygon(belling_outline, XG) | in_polygon(amundsen_outline, XG) );
             temp=read_gmt_vectors('/Volumes/ice1/ben/Cryosat/clip_files/coastal_amundsen.gmt');
             coastal=ll2ps(temp.ll(:,2), temp.ll(:,1));
             XG=XG(in_polygon(coastal,XG));
             disp('RUNNING PIG ONLY:'); XG=1e3*(-1600-250i)
             region_params.E_RMS_d2z0dx2=2e-6;
             % changed from 6e-8 to 1e-8, rerun
             region_params.E_RMS_d3zdx2dt=1e-8;
             region_params.E_RMS_d2zdt2=0.5;
         end
         if 0
             % near the big lake, testing how much to relax the parameters
             % using the roughness map.
             XG=(-1450+1i*-450)*1e3;
             region_params.E_RMS_d2z0dx2=2e-7;  % marginally strong enough to suppress small peaks
             region_params.E_RMS_d3zdx2dt=6e-8;
             region_params.E_RMS_d2zdt2=1;
             region_params.power_fit=false;
            % For testing:
             %region_params.skip_MOA_curv= true;
             
         end
         
        if 0
            % THW GL
            XG= -1.5500e+06 - 4.5000e+05i + 5e4;
            region_params.E_RMS_d2z0dx2=2e-7/10;
            region_params.E_RMS_d3zdx2dt=6e-8/10;
            region_params.E_RMS_d2zdt2=0.2;
        end
        
        if 0
            % THW GL, relaxed for real
            XG= -1.5500e+06 - 4.5000e+05i + 4e4;
            region_params.E_RMS_d2z0dx2=2e-7;
            region_params.E_RMS_d3zdx2dt=6e-8;
            region_params.E_RMS_d2zdt2=0.2;
        end
        
        
        if 0
            % the big lake
            XG=-1.4110e+06+1i*-4.3033e+05;  
            %XG=(-1425-1i*400)*1e3;
            region_params.E_RMS_d2z0dx2=2e-7;
            region_params.E_RMS_d3zdx2dt=6e-8;
            region_params.E_RMS_d2zdt2=1;
        end
        if 0
            %XG=(-1400+1i*-425)*1e3 +
            XG=(-1450+1i*-375+[-25 0 25])*1e3;
            region_params.E_RMS_d2z0dx2=2e-7;
            region_params.E_RMS_d3zdx2dt=6e-8;
            region_params.E_RMS_d2zdt2=1;
            region_params.est_errors=false;
        end
        
end
XG=unique(XG);


%--------------------------------------------------------------------------
function [S, D]=fit_AA(sub_region, params_in)

dT=0.25;
dx0=500; dx=1e3;

%delta=2.5e4; W=delta+2e4;
delta=5e4; W=delta+2e4;

baseline_b=exist('params_in','var') && isfield(params_in,'baseline_b') && params_in.baseline_b;
if ~baseline_b
    baseline_c=true;
    params_in.baseline_c=true;
end

CS_data_new='/Volumes/insar6/ben/Cryosat/';

%try the new Helm picks
data_file_POCA=[CS_data_new,'POCA_h5_C/south/AA_Helm_v2c_sigma4/master_index_h5.mat'];
data_file_sw=[CS_data_new,'SW_h5_C/south/AA_Helm_v2c_sigma4/master_index_h5.mat'];

%out_dir='/Volumes/ice1/ben/Cryosat/fits/Amundsen_combo_v2/';
%out_dir='/Volumes/ice1/ben/Cryosat/fits/Amundsen_lakes_SwathBiasOnly_relaxed/';
%out_dir='/Volumes/ice1/ben/Cryosat/fits_C/Amundsen_lakes_SwathBiasOnly_relaxed_v2/';
out_dir='/Volumes/ice1/ben/Cryosat/fits_C/Helm_picks_relazed_z0/';
if strcmp(sub_region,'test')
    % HELM
    %out_dir='/Volumes/ice1/ben/Cryosat/fits_C/testing_helm/';
    out_dir='/Volumes/ice1/ben/Cryosat/fits_C/testing_helm_data_Zmodel_fit/';
end
     
%out_dir='/Volumes/ice1/ben/Cryosat/fits/Amundsen_lakes_SwathBiasOnly_relaxed_coarse';
%out_dir='/Volumes/ice1/ben/Cryosat/fits/Belling_POCA_July30/';
%out_dir='/Volumes/ice1/ben/Cryosat/fits/test_z0_2e-7_dz_6e-8/';
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
% OLD:
%DEM_file='/Volumes/insar5/gmap/OIB_data/AA/fits/100_km/3e-5_dz_SScorr_firncorr_test_nov16_stereo_GC_50cm_error_BiasCorr/dZ_grid.mat';
%Helm:
DEM_file='/Volumes/ice1/ben/helm_DEM/AA/helm_DEM_with_EGM2008.tif';
% Newer:
%DEM_file='/Volumes/insar5/gmap/OIB_data/AA/fits/100_km/3e-5_dz_SScorr_firncorr_Nov2015_stereo_r634_50cm_error_DEMSlope_1m_100km_v2_initial/dZ_grid.mat';
%Zmodel?
%DEM_file='/Volumes/ice1/ben/Cryosat/fits_C/Zmodel_Dec_1_2016.mat';

[XG, region_params]=define_region(sub_region, delta, data_file_POCA );
[~, ind]=sort(abs(XG-mean(XG))); XG=XG(ind);
YG=imag(XG); XG=real(XG);

if exist('params_in','var') && isfield(params_in,'X0')
    XG=real(params_in.X0);
    YG=imag(params_in.X0);
    OverWrite=true;
else
    OverWrite=false;
end
params=struct('W', W,'dx', dx, 'dx0', dx0, 'dt', dT, 'TR', [10.5 16.5], 'swath', true,...
    'POCA', true, 'DEM_file', DEM_file,'est_errors', false,'time_zero_season', 5, 'out_dir', out_dir);%, 'calc_PS_bias', true);
%params for lakes: probably lax for the d2zdt2 parameter.  Set these as the
%default
params.E_RMS_d2z0dx2=2e-7;
params.E_RMS_d3zdx2dt=6e-8;
params.E_RMS_d2zdt2=1;

% map region params to params, overwrite defaults
if isstruct(region_params);
    f=fieldnames(region_params);
    for kF=1:length(f);
        params.(f{kF})=region_params.(f{kF});
    end
end

% define sigma structure:
line_spacing_z0=1e3;
line_spacing_dz=1e3;
sigma=struct( 'smooth_dz',params.E_RMS_d3zdx2dt*W, 'smooth_z0', params.E_RMS_d2z0dx2*W, 'same_season_dz', 2*W,...
    'smooth_season_dz', params.E_RMS_d2zdt2*W,'smooth_bias', .001*W, 'zero_bias', .001*W, 'line_spacing_z0', line_spacing_z0, 'line_spacing_dz', line_spacing_dz );
%sigma.time_gap=0.5;
sigma.time_gap=2;
%sigma.smooth_bias=sigma.smooth_bias/1e8;
sigma.zero_bias=2*W;
sigma.smooth_bias=0.01*params.E_RMS_d2z0dx2*W;
sigma.seasonal_cycle=0.01;
params.sigma=sigma;

params.calc_roll_bias=false;%true;
params.calc_PS_bias=true;
params.calc_perswath_bias=true;
if isempty(data_file_sw); 
    params.calc_roll_bias=false;%true;
    params.calc_PS_bias=false;
    params.calc_perswath_bias=false;
    params.calc_AD_bias=true;
end

if params.est_errors==0; disp('skipping errors!!!');end
 
load /Volumes/ice1/ben/Cryosat/Zmodel_Dec_10_2014
params.DEM=struct('x',Zmodel.x, 'y', Zmodel.y, 'z', Zmodel.z0);
params.dZ_model=struct('x', Zmodel.x, 'y', Zmodel.y, 'z', Zmodel.z, 't0', Zmodel.t0);

for k=1:length(XG(:))
    save_name=sprintf('%s/fit_%d_%d.mat', params.out_dir, round(XG(k)/1e3), round(YG(k)/1e3));
    if (~exist(save_name,'file') && lockfile_tool('check', save_name)==0) || OverWrite
        lockfile_tool('lock', save_name);
        fprintf(1,'---------working on %s----------\n', save_name);
        %com.mathworks.mde.desk.MLDesktop.getInstance.getMainFrame.setTitle(['Matlab: proc: ', save_name]);
        D=load_both(XG(k)+1i*YG(k), {data_file_POCA, data_file_sw}, params);
        
        % if the mask exists, interpolate it to the data points
        if isfield(params,'mask_file') && ~isempty(params.mask_file);
            mask=read_geotif_xy(params.mask_file, real(params.X0)+[-1.1 1.1]*params.W/2, imag(params.X0)+[-1.1 1.1]*params.W/2);
            mask.z=dilate_mask(mask.z~=0, 3);
            D.mask_at_pts=interp2(mask.x, mask.y, double(mask.z~=0), D.x, D.y,'*nearest')~=0;
        end
        
        D.year=(D.time-datenum('jan 1 2000'))/365;
        
        M=struct('dx', params.dx, 'dy', params.dx, 'dt', params.dt, 'XR', XG(k)+[-1 1]*W/2,'YR', YG(k)+[-1 1]*W/2, 'TR', params.TR,'dx0', params.dx0,'dy0', params.dx0,'time_zero_season', params.time_zero_season, ...
            'sigma', params.sigma);
        
        if isempty(D.x);
            S=[];
           continue;
        end                
        [S, D]=fit_section(D, M, params);
        save(save_name, 'S','D','M','params');
        lockfile_tool('unlock', save_name);
        fprintf(1,'---------done with %s----------\n\n\n', save_name);
    end
end

%------------------------------------------------
function [S, D]=fit_GL(sub_region, params_in)

dT=0.25;
dx0=250; dx=1e3;
delta=5e4; W=delta+2e4;
DEM_file='/Volumes/insar9/ben/GIMP2_DEM/v2fit_120.tif';

CS_data_new='/Volumes/insar6/ben/Cryosat/';
data_file_POCA=[CS_data_new,'POCA_h5_C/north/gimp_v1c_sigma4/master_index_h5.mat'];
data_file_sw=[CS_data_new,'SW_h5_C/north/gimp_v1c_sigma4/master_index_h5.mat'];

out_dir='/Volumes/ice1/ben/Cryosat/fits_C/gimp2_picks_relazed_z0/';
if strcmp(sub_region,'test')
     out_dir='/Volumes/ice1/ben/Cryosat/fits_C/testing_gimp_picks';
end
   
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
 
[XG, region_params]=define_region(sub_region, delta, data_file_POCA );
[~, ind]=sort(abs(XG-mean(XG))); XG=XG(ind);
YG=imag(XG); XG=real(XG);

if exist('params_in','var') && isfield(params_in,'X0')
    XG=real(params_in.X0);
    YG=imag(params_in.X0);
    OverWrite=true;
else
    OverWrite=false;
end
params=struct('W', W,'dx', dx, 'dx0', dx0, 'dt', dT, 'TR', [10.5 16.5], 'swath', true,...
    'POCA', true, 'DEM_file', DEM_file,'est_errors', false,'time_zero_season', 5, 'out_dir', out_dir);%, 'calc_PS_bias', true);
%params for lakes: probably lax for the d2zdt2 parameter.  Set these as the
%default
params.E_RMS_d2z0dx2=2e-7;
params.E_RMS_d3zdx2dt=6e-8;
params.E_RMS_d2zdt2=1;

% map region params to params, overwrite defaults
if isstruct(region_params);
    f=fieldnames(region_params);
    for kF=1:length(f);
        params.(f{kF})=region_params.(f{kF});
    end
end

% define sigma structure:
line_spacing_z0=1e3;
line_spacing_dz=1e3;
sigma=struct( 'smooth_dz',params.E_RMS_d3zdx2dt*W, 'smooth_z0', params.E_RMS_d2z0dx2*W, 'same_season_dz', 2*W,...
    'smooth_season_dz', params.E_RMS_d2zdt2*W,'smooth_bias', .001*W, 'zero_bias', .001*W, 'line_spacing_z0', line_spacing_z0, 'line_spacing_dz', line_spacing_dz );
%sigma.time_gap=0.5;
sigma.time_gap=2;
%sigma.smooth_bias=sigma.smooth_bias/1e8;
sigma.zero_bias=2*W;
sigma.smooth_bias=0.01*params.E_RMS_d2z0dx2*W;
sigma.seasonal_cycle=0.01;
params.sigma=sigma;

params.calc_roll_bias=false;%true;
params.calc_PS_bias=true;
params.calc_perswath_bias=true;
 

if params.est_errors==0; disp('skipping errors!!!');end

for k=1:length(XG(:))
    save_name=sprintf('%s/fit_%d_%d.mat', params.out_dir, round(XG(k)/1e3), round(YG(k)/1e3));
    if (~exist(save_name,'file') && lockfile_tool('check', save_name)==0) || OverWrite
        lockfile_tool('lock', save_name);
        fprintf(1,'---------working on %s----------\n', save_name);
        D=load_both(XG(k)+1i*YG(k), {data_file_POCA, data_file_sw}, params);
        
        % if the mask exists, interpolate it to the data points
        if isfield(params,'mask_file') && ~isempty(params.mask_file);
            mask=read_geotif_xy(params.mask_file, real(params.X0)+[-1.1 1.1]*params.W/2, imag(params.X0)+[-1.1 1.1]*params.W/2);
            mask.z=dilate_mask(mask.z~=0, 3);
            D.mask_at_pts=interp2(mask.x, mask.y, double(mask.z~=0), D.x, D.y,'*nearest')~=0;
        end
        
        % note: changed from 365 to 365.25 3/1/2018
        D.year=(D.time-datenum('jan 1 2000'))/365.25;
        
        M=struct('dx', params.dx, 'dy', params.dx, 'dt', params.dt, 'XR', XG(k)+[-1 1]*W/2,'YR', YG(k)+[-1 1]*W/2, 'TR', params.TR,'dx0', params.dx0,'dy0', params.dx0,'time_zero_season', params.time_zero_season, ...
            'sigma', params.sigma);
        
        if isempty(D.x);
            S=[];
           continue;
        end                
        [S, D]=fit_section(D, M, params);
        save(save_name, 'S','D','M','params');
        lockfile_tool('unlock', save_name);
        fprintf(1,'---------done with %s----------\n\n\n', save_name);
    end
end




%-------------------------------------------------
function D=load_swath_data(X0, W, data_file)

load(data_file)
INDEX=master_index;

dx=INDEX.dx;
xg=(floor((real(X0)-W/2)/dx):ceil((real(X0)+W/2)/dx))*dx;
yg=(floor((imag(X0)-W/2)/dx):ceil((imag(X0)+W/2)/dx))*dx;
[xg, yg]=meshgrid(xg, yg);

D=index_point_data_h5('read_from_index',   xg(:)+1i*yg(:),  INDEX, {'x','y','time','h', 'power','coherence','AD','error_composite', 'R_POCA', 'ambiguity','burst','abs_orbit', 'block_h_spread','count','phase','R_offnadir','range_surf','seg_ind','dRange_POCA'}, []);
if any(D.count>1)
    D=index_struct(D, D.power > 1e-17 & D.power < 1e-13 & D.error_composite==0 & D.count > 3 & D.block_h_spread < 15);
else
    D=index_struct(D, D.power > 1e-17 & D.power < 1e-13 & D.error_composite==0);
end

if 0
    D=index_struct(D, ~ismember(D.abs_orbit, [ 1613 16900 16907 18117 18122, 18123, 22468, 24182 24183, 24185, 24190 24191 24192 24177 26305 26306]));
    
    load /Volumes/ice1/ben/Cryosat/BadOrbits_5m_April_4_2016.mat
    D=index_struct(D, ~ismember(D.abs_orbit, BadOrbits));
end
% t_bins=min(D.time):(365/10):max(D.time)+(365/10);
% for k=1:length(t_bins)-1;
%     temp{k}=find(D.time >=t_bins(k) & D.time < t_bins(k+1));
%     if ~isempty(temp{k});
%         temp{k}=temp{k}(blockmedian(D.x(temp{k})+1i*D.y(temp{k}), D.h(temp{k}), 250));
%     end;
% end
% D=index_struct(D, cat(1, temp{:}));

%--------------------------------------------------------------------------
function D=load_POCA_data(X0, W, data_file)

if strcmp(data_file(end-2:end),'mat')
    load(data_file)
    INDEX=master_index;
else
    INDEX=index_point_data_h5('read_index', data_file);
end
dx=INDEX.dx;
xg=(floor((real(X0)-W/2)/dx):ceil((real(X0)+W/2)/dx))*dx;
yg=(floor((imag(X0)-W/2)/dx):ceil((imag(X0)+W/2)/dx))*dx;
[xg, yg]=meshgrid(xg, yg);
%D=index_point_data_h5('read_h5_file', data_file,  xg(:)+1i*yg(:),  {'x','y','time','surf_height', 'backsc_sig','AD','x_sc','y_sc', 'h_sc','Qflag'}, []);

if strcmp(data_file(end-2:end),'mat')
    D=index_point_data_h5('read_from_index',   xg(:)+1i*yg(:),  INDEX, {'x','y','time','h', 'power','coherence','AD','error_composite', 'ambiguity','burst','abs_orbit','phase','range_surf'}, []);
else
    D=index_point_data_h5('read_h5_file', data_file,  xg(:)+1i*yg(:),  {'x','y','time','h', 'power','coherence','AD','error_composite', 'ambiguity','pulse_num','abs_orbit','phase','range_surf'}, []);
    D.burst=D.pulse_num;
    D=rmfield(D,'pulse_num');
end

%--------------------------------------------------------------------------
function D=filter_POCA_data(D)

% baseline B:
%els=D.power> 5e-15 & D.coherence> 800 & D.error_composite==0 & D.power < 1e-12;
% Baseline C:
els=D.power > 1e-16 & D.error_composite == 0 & D.power < 1e-12;

els(abs(D.time-datenum('28 jul 2010 18:05'))<10/60/24)=false;

els=els & ~ismember(D.abs_orbit, [ 4069 16900:16907 18117:18123, 22468, 24182:24185, 24190:24192 24177 26306]  );

load /Volumes/ice1/ben/Cryosat/BadOrbits_5m_April_4_2016.mat
els=els & ~ismember(D.abs_orbit, BadOrbits);

if isfield(D,'mask_at_pts');
    els=els&D.mask_at_pts;
end
D=index_struct(D, els);

%------------------------------------------------------------------------
function D=load_both(X0, data_files, params)

% load the POCA data
D{1}=load_POCA_data(X0, params.W, data_files{1});

% assign errors:
D{1}.sigma=ones(size(D{1}.x));
D{1}=filter_POCA_data(D{1});
D{1}.swath=false(size(D{1}.time));

try 
    temp=load(params.DEM_file);
catch
    temp.dZ=read_geotif_xy(params.DEM_file, real(X0)+[-0.6 0.6]*params.W, imag(X0)+[-0.6 0.6]*params.W);
    temp.dZ.z0=temp.dZ.z;
    temp.dZ=rmfield(temp.dZ,'z');
end
[gx, gy]=gradient(temp.dZ.z0, temp.dZ.x, temp.dZ.y);
DEM_slope_mag=abs(interp2(temp.dZ.x, temp.dZ.y, gx+1i*gy, D{1}.x, D{1}.y));
D{1}.sigma=max(0.5, 50*DEM_slope_mag+ max(0, -0.64*(log10(D{1}.power)+14)));


% load the swath data
if ~isempty(data_files{2});
    D{2}=load_swath_data(X0, params.W, data_files{2});
    
    %D{2}.sigma=max(1, D{2}.block_h_spread);
    % use parameters from twoD_PS_error_analysis_swath
    D{2}.sigma=min(5, max(1, 0.95 -.4816*(log10(D{2}.power)+14) + 1.12*sqrt(D{2}.block_h_spread)));
    
    D{2}.swath=true(size(D{2}.time));
    D{2}=index_struct(D{2},D{2}.R_POCA >1000);
    
    % jam the data together
    D1_count=length(D{1}.time);
    D2_count=length(D{2}.time);
    f2=fieldnames(D{2});
    for kF=1:length(f2);
        if isfield(D{1}, f2{kF});
            D{1}.(f2{kF})=[D{1}.(f2{kF}); D{2}.(f2{kF})];
        else
            D{1}.(f2{kF})=[  NaN(D1_count, 1); D{2}.(f2{kF})];
        end
    end
    
    f1=fieldnames(D{1});
    for kF=1:length(f1);
        if ~isfield(D{2}, f2{kF});
            D{1}.(f1{kF})=[D{1}.(f2{kF}); NaN(D2_count,1)];
        end
    end
end
D=D{1};

%----------------------------------------------
function eq_key=define_equation_types
% local version of est_dh_dt_fd's function 'define_equation_types', adding
% types for additional equations

eq_key=est_dh_dt_fd('define_equation_types');
eq_key.roll_fit=31;
eq_key.power_fit=21;
eq_key.perswath_fit=32;

%--------------------------------------------------------------------------
function [S, D]=fit_section(D, M, params)

% given data structure 'D'  model structure 'M' and param structure
% 'params'
% set up the fitting matrices and iterate
DOPLOT=false; %true;
info=default_info;

% define equation typess
eq_key=define_equation_types;

%assemble the matrices for the fit, not including the biases
[grids, G, Gc, Cvals, TOC, D]=est_dh_dt_fd('prep_fit', D, M);
if isempty(D.x);
    S=[]; 
    D=[]; 
    return
end

% read in the MOA curvature, add it in quadrature to the estimated
% curvature scale
%Mcurv=read_geotif_xy('/Volumes/ice1/ben/MOA/MOA_curve_burned.tif', range(grids.z0.ctrs{2})+[-1e4 1e4], range(grids.z0.ctrs{1})+[-1e4 1e4]);
if ~isfield(params,'skip_MOA_curv') && info.hemisphere==-1;
    Mcurv=read_geotif_xy('/Volumes/ice1/ben/MOA/MOA_curv.tif', range(grids.z0.ctrs{2})+[-1e4 1e4], range(grids.z0.ctrs{1})+[-1e4 1e4]);
    these=TOC.eqn_type==eq_key.smooth_z0;
    temp=Gc(these,TOC.cols.z0)~=0;
    [gy, gx]=ndgrid(grids.z0.ctrs{1:2});
    xloc=temp*gx(:); xloc=xloc./sum(temp,2);
    yloc=temp*gy(:); yloc=yloc./sum(temp,2);
    Curve_est=double(interp2(Mcurv.x, Mcurv.y, Mcurv.z, xloc, yloc));
    Curve_est(~isfinite(Curve_est))=0;
    Cvals.Gc(these)=Cvals.Gc(these)+(Curve_est*params.W*10).^2;
    
    % do the same for the flat_z0 params
    these=TOC.eqn_type==eq_key.flat_z0;
    temp=Gc(these,TOC.cols.z0)~=0;
    xloc=temp*gx(:); xloc=xloc./sum(temp,2);
    yloc=temp*gy(:); yloc=yloc./sum(temp,2);
    Curve_est=double(interp2(Mcurv.x, Mcurv.y, Mcurv.z, xloc, yloc));
    Curve_est(~isfinite(Curve_est))=0;
    Cvals.Gc(these)=Cvals.Gc(these)+(Curve_est*params.W*1000*10).^2;
end
[count, trange, mask_bin]=count_CS_data(D, grids,'all');
mask=count>5 & trange >2;
mask=dilate_mask(~dilate_mask(~mask, 3), 4);
if mean(mask(:))<0.1;
    disp('not enough data to run a fit (< 10% of grid cells); returning');
    S=[]; D=[]; return;
end

bias_count=0;

% add the equations for the AD bias
if isfield(params,'calc_AD_bias') && params.calc_AD_bias
    [G_fit_bias, Gc_bias, bias_grids, TOC_bias]=est_dh_dt_fd('bias_interp_mtx', M, D, -0.5*(D.AD==-1) + 0.5*(D.AD==1), eq_key);
    bias_count=bias_count+1;
    bias_eq(bias_count)=struct('name', 'AD_bias','G', G_fit_bias,'Gc', Gc_bias, 'TOC', TOC_bias);
elseif isfield(params,'calc_PS_bias') && params.calc_PS_bias
    %    [G_fit_bias, Gc_bias, bias_grids, TOC_bias]=est_dh_dt_fd('bias_interp_mtx', M, D, -0.5*(D.swath==0) + 0.5*(D.swath==1), eq_key);
    [G_fit_bias, Gc_bias, bias_grids, TOC_bias]=est_dh_dt_fd('bias_interp_mtx', M, D, double(D.swath==1), eq_key);
    bias_count=bias_count+1;
    bias_eq(bias_count)=struct('name', 'PS_bias','G', G_fit_bias,'Gc', Gc_bias, 'TOC', TOC_bias);
end
if exist('bias_grids','var');
    grids.bias=bias_grids.bias;
end

% add the season bias
if isfield(params.sigma, 'seasonal_cycle');
    bias_count=bias_count+1;
    TOC_season=struct('eqn_id', [0; 0],'eqn_type', [1; 1]*eq_key.seasonal_cycle, 'node_num', [0; 0] );
    bias_eq(bias_count)=struct('name', 'seasonal_cycle','G', [sin(2*pi*D.year) cos(2*pi*D.year)],'Gc', eye(2), 'TOC', TOC_season);
end

% add a roll-bias estimate for each orbit
if isfield(params,'calc_roll_bias') && params.calc_roll_bias
    dYdPhi=7.5e5*.022084/2/pi/1.17;
    [uOrb, ~, iOrb]=unique(D.abs_orbit(D.swath==1).*D.AD(D.swath==1));
    G_roll=sparse(find(D.swath==1), iOrb, double(D.phase(D.swath==1)+2*pi*D.ambiguity(D.swath==1))*dYdPhi, length(D.x), length(uOrb));
    Gc_roll=speye(size(G_roll,2));
     
    TOC_roll_bias=struct('eqn_id', eq_key.roll_fit*ones(size(G_roll,2),1), ...
        'eqn_type', eq_key.roll_fit*ones(size(G_roll,2),1),...
        'node_num', zeros(size(G_roll,1),1));
         
    bias_count=bias_count+1;
    bias_eq(bias_count)=struct('name', 'roll_fit',...
        'G', G_roll,'Gc', Gc_roll, 'TOC', TOC_roll_bias);
end

% add a uniform bias estimate for each orbit
if isfield(params,'calc_perswath_bias') && params.calc_perswath_bias
    [uOrb, ~, iOrb]=unique(D.abs_orbit(D.swath==1).*D.AD(D.swath==1));
    G_perswath=sparse(find(D.swath==1), iOrb, ones(sum(D.swath==1),1), length(D.x), length(uOrb));
    Gc_perswath=speye(size(G_perswath,2));
     
    TOC_perswath=struct('eqn_id', eq_key.perswath_fit*ones(size(G_perswath,2),1), ...
        'eqn_type', eq_key.perswath_fit*ones(size(G_perswath,2),1),...
        'node_num', zeros(size(G_perswath,1),1));
         
    bias_count=bias_count+1;
    bias_eq(bias_count)=struct('name', 'perswath_fit',...
        'G', G_perswath,'Gc', Gc_perswath, 'TOC', TOC_perswath);
end

% add another (unconstrained) parameter for power
if isfield(params,'power_fit') && params.power_fit
    TOC_power_fit=struct('eqn_type', eq_key.power_fit,'eqn_id', eq_key.power_fit,'node_num', zeros(0, 1));
    bias_count=bias_count+1;
    bias_eq(bias_count)=struct('name', 'power_fit',...
        'G', double(log10(D.power)-median(log10(D.power))),'Gc', zeros(0,1), 'TOC', TOC_power_fit);
end
    
    
% add the bias equations to G, Gc, and TOC
last_bias_col=size(G,2);
for k=1:length(bias_eq);
    these_cols=last_bias_col+(1:size(bias_eq(k).G,2));
    G=[G, bias_eq(k).G];
    Gc=[Gc sparse( size(Gc,1), size(bias_eq(k).G,2)); ...
        sparse( size(bias_eq(k).Gc,1), size(Gc,2)), bias_eq(k).Gc];
    TOC.cols.(bias_eq(k).name)=these_cols;
    TOC=cat_struct(TOC, bias_eq(k).TOC, 1, {'eqn_id','eqn_type','node_num'});
    last_bias_col=these_cols(end);
end

% add Cvals for the biases to the Cvals for the data
Cvals.Gc(TOC.eqn_type==eq_key.smooth_bias)=M.sigma.smooth_bias^2;
Cvals.Gc(TOC.eqn_type==eq_key.flat_bias)=M.sigma.smooth_bias.^2*M.sigma.line_spacing_z0.^2;
Cvals.Gc(TOC.eqn_type==eq_key.zero_bias)=M.sigma.zero_bias.^2;
Cvals.Gc(TOC.eqn_type==eq_key.seasonal_cycle)=M.sigma.seasonal_cycle.^2;
if isfield(params,'calc_roll_bias') && params.calc_roll_bias
    Cvals.Gc(TOC.eqn_type==eq_key.roll_fit)=(1/12000).^2;
end
if isfield(params,'calc_perswath_bias') && params.calc_perswath_bias
    Cvals.Gc(TOC.eqn_type==eq_key.perswath_fit)=10^2;
end

d_c=zeros(size(Gc,1),1);

% iterate the fit
good=true(size(G,1),1) & mask(mask_bin);
if ~any(good);
    S=[]; D=[];
    return
end
curve_edit=good;
N_edit_iterations=15;
m=zeros(size(G,2),1);
N_clip=3;
count_sw=count_CS_data(index_struct(D,   D.swath==1), grids,'all');
count_POCA=count_CS_data(index_struct(D, D.swath==0), grids,'all');
for k=1:N_edit_iterations;
    % assemble the C matrix
    
    if k>1;
        C_extra=1./((1-(rs(good)./threshold).^2)).^2;
    else
        C_extra=ones(sum(good),1);
    end
    
    C=spdiags([Cvals.G(good,:).*C_extra;  Cvals.Gc], 0, sum(good)+size(Gc,1), sum(good)+size(Gc,1));
    m_last=m;
    % run the fit
    tic;
    m=my_lscov([G(good,:); Gc], double([D.h(good); d_c]), C);this_dt=toc;
    r=D.h-G*m;
    rs=r./sqrt(Cvals.G);
    [dR_med, local_thresh]=local_robust_stats(D.x+1i*D.y, double(r), good,  grids.z0.ctrs{2}(1)+1i*grids.z0.ctrs{1}(1), 5e3);
    %rs=dR_med./sqrt(Cvals.G);
    sigma_hat=iqr(rs(good))/2;   
    last_good=good;
     if k==1;
         threshold=max(20, N_clip*sigma_hat);
     else
         threshold=max(N_clip, N_clip*sigma_hat);
     end
     %good=mask(mask_bin) & abs(rs) < threshold*max(1, local_thresh);

     last_good=good;
     good=abs(rs)<threshold & mask(mask_bin);
%      if k>5
%          good=good & abs(r) < 15;
%      end
    R2=est_dh_dt_fd('parse_residuals', Gc*m, Cvals.Gc, TOC);
    R2.data=sum(rs(last_good).^2);
    
    if k==1;
        fprintf(1, 'it=%d, sigma_hat=%4.2f, threshold=%4.2f, P_g=%3.2f, dt=%2.1f\n', k, sigma_hat, threshold,  mean(good)*100, this_dt);
    else
        fprintf(1, 'it=%d, sigma_hat=%4.2f, threshold=%4.2f, P_g=%3.2f%%, dt=%2.1f,  max(dm_dz)=%f, max(dm_z0)=%f\n', ...
            k, sigma_hat, threshold,  mean(good)*100, this_dt,  max(abs(m_last(TOC.cols.dz)-m(TOC.cols.dz))),max(abs(m_last(TOC.cols.z0)-m(TOC.cols.z0)))  );
    end
    %fprintf('  R2c_z0=%f, R2c_dz=%f\n', R2.z0, R2.dz);
    %fprintf('  R2c_total=%6.1f, R2_data=%6.1f, R2_data/R2_total=%4.2f, R2_data_scaled/R2c_total=%f\n', R2.total, R2.data, R2.data/(R2.total+R2.data), R2.data/(R2.total.*params.CEQ_wt.^2));
    
    if DOPLOT;
        dSi=est_dh_dt_fd('parse_fit', G, Gc, m_last-m, D.h, TOC, grids);
        Si=est_dh_dt_fd('parse_fit', G, Gc, m, D.h, TOC, grids);
        figure(k); clf; set(gcf,'position', [113    77   590   944]);
        hh=cheek_by_jowl(4, 2, [0 0 1 1 ])'; 
        axes(hh(1));  hs=surf(grids.z0.ctrs{2}, grids.z0.ctrs{1}, Si.z0); shading interp; set(gca,'dataaspectratio', [ 1 1 .01]);  hl=light; view([0 90]); material dull; axis tight;
        %axes(hh(3)); imagesc(std(Si.dz, [], 3)); axis xy equal; colorbar;
        axes(hh(2)); imagesc((Si.dz(:,:,end-1)-Si.dz(:,:,2))/diff(grids.dz.ctrs{3}([2 end]))); caxis([-2 2]); axis xy equal; colorbar;
        last_count_sw=count_sw;
        last_count_POCA=count_POCA;
        count_sw=count_CS_data(index_struct(D, good & D.swath==1), grids,'all');
        count_POCA=count_CS_data(index_struct(D, good & D.swath==0), grids,'all');
        axes(hh(3)); imagesc(dSi.z0); axis xy equal; caxis([-5 5]); 
        axes(hh(4)); imagesc(std(dSi.dz, [], 3)); axis xy equal;
        axes(hh(5)); imagesc(count_sw); caxis([1 100]); axis xy equal tight; 
        axes(hh(6)); imagesc(count_POCA); caxis([1 10]); axis xy equal tight;
        axes(hh(7)); imagesc(count_sw - last_count_sw); caxis([-10 10]); axis xy equal tight; 
        axes(hh(8)); imagesc(count_POCA - last_count_POCA); caxis([-10 10]); axis xy equal tight; 
        colormap([0 0 0; jet(127)*.7+.3]);
        drawnow
    end

    z0=reshape(m(TOC.cols.z0), grids.z0.dims);
    temp=sqrt(conv2(z0, [-1 2 -1]/2,'same').^2+conv2(z0, [-1 2 -1]'/2,'same').^2);
    temp(1,:)=0; temp(end,:)=0; temp(:,1)=0; temp(:,end)=0;
    [rr,cc]=find(temp>20);
    curve_edit=true(size(good));
    for kSpike=1:length(rr);
        xy0=[grids.z0.ctrs{2}(cc(kSpike))+1i*grids.z0.ctrs{1}(rr(kSpike))];
        close_x= abs(D.x+1i*D.y-xy0)<3*diff(grids.z0.ctrs{2}(1:2));
        if any(good & close_x);
            good_temp=abs(r)<0.8*max(abs(r(good & close_x)));
            curve_edit(close_x & ~good_temp)=false;
        end
    end
    fprintf(1, 'curvature edit removes %d otherwise good values, %d total \n', sum(good &~curve_edit), sum(~curve_edit));
    
    good(~curve_edit)=0;
    
    if all(last_good==good) || max(abs(m_last(TOC.cols.dz)-m(TOC.cols.dz)))<1e-1;
        fprintf(1, 'done iterating- acceptance criteria not changing or model converged\n');
        good=last_good;  % rewind- we're
        break
    end
    
    % give the ierations a chance to converge, then try to remove the worst swaths
    if k>=4 && isfield(params,'calc_perswath_bias') && params.calc_perswath_bias
        BadPerSwathOrbit=uOrb(abs(m(TOC.cols.perswath_fit))>10);
    end
    if exist('BadPerSwathOrbit','var')
        good=good & ~ismember(D.abs_orbit.*D.AD, BadPerSwathOrbit);
    end
    
    if false
        % flag swath segments that have RMS(rs)> 5
        if k>=4
            uOBS=unique([D.abs_orbit(D.good & D.swath), D.burst(D.good& D.swath), D.seg_ind(D.good& D.swath)],'rows');
            [~, ii]=ismember([D.abs_orbit, D.burst, D.seg_ind, D.swath], [uOBS, ones(size(uOBS,1),1)],'rows');
            [~, jj]=ismember(ii, 1:length(ii));
            DS_mat=sparse(jj(ii~=0), ii(ii~=0), ones(sum(ii~=0),1), length(uOBS), length(D.x));
            SS_seg=DS_mat*(D.good.*rs).^2;
            RMS_seg=SS_seg;
            N_seg=sum(DS_mat,2);
            RMS_seg(N_seg>0)=RMS_seg(N_seg>0)./N_seg(N_seg>0);
            RMS_seg=sqrt(RMS_seg);
        end
    end
    
end

% report the results
S=est_dh_dt_fd('parse_fit', G, Gc, m, D.h, TOC, grids);
if isfield(TOC.cols,'PS_bias');
    %S.bias=reshape(m(TOC.cols.AD_bias), grids.bias.dims);
    S.bias=reshape(m(TOC.cols.PS_bias), grids.bias.dims);
end

if isfield(TOC.cols,'AD_bias');
     S.bias=reshape(m(TOC.cols.AD_bias), grids.bias.dims);
end
   
S.grids=grids;
S.seasonal_cycle=m(TOC.cols.seasonal_cycle);

[S.count, S.t_range, ~, S.count_z0]=count_CS_data(index_struct(D, good), grids,'all');
[S.count_by_season]=count_CS_data(index_struct(D, good), grids,'by_season');

if isfield(TOC.cols, 'power_fit')
    S.power_fit=m(TOC.cols.power_fit);
end

D.good=good;
D.d_est=G*m;
if isfield(params,'est_errors') && params.est_errors
    C_extra=1./((1-(rs(good)./threshold).^2)).^2;
    C=spdiags([Cvals.G(good,:).*C_extra;  Cvals.Gc], 0, sum(good)+size(Gc,1), sum(good)+size(Gc,1));
    
    tic; [~, sigma_m, mse]=my_lscov([G(good,:); Gc], double([D.h(good); d_c]), C);toc
    sigma_m=sigma_m*sqrt(1/mse);
    Se=est_dh_dt_fd('parse_fit', G, Gc, sigma_m, D.h, TOC, grids);
    S.sigma_dz=Se.dz;
    S.sigma_z0=Se.z0;
end

if isfield(params,'calc_roll_bias') && params.calc_roll_bias
    S.roll_bias.AbsOrb=uOrb;
    S.roll_bias.bias=m(TOC.cols.roll_fit);
    if isfield(params,'est_errors') && params.est_errors
        S.roll_bias.sigma=sigma_m(TOC.cols.roll_fit);
    end
end

if isfield(params,'calc_perswath_bias') && params.calc_perswath_bias
    S.perswath_bias.AbsOrb=uOrb;
    S.perswath_bias.bias=m(TOC.cols.perswath_fit);
    if isfield(params,'est_errors') && params.est_errors
        S.perswath_bias.sigma=sigma_m(TOC.cols.perswath_fit);
    end
end

m1=zeros(size(m));
m1(TOC.cols.dz)=m(TOC.cols.dz);
D.dz=G*m1;

m1=zeros(size(m));
m1(TOC.cols.z0)=m(TOC.cols.z0);
D.z0=G*m1;









