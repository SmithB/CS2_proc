function D_out=locate_CS_returns_BD(D_in, D, CP, hemisphere, params)

burst_num=D_in.burst;
sfc_bin=D_in.samp;

 
sfc_phase=D_in.phase;

[RefEll.a,RefEll.b,RefEll.e2,RefEll.finv]=refell('WGS84');

% setup the polar_stereographic projection to match the hemisphere 
if ~exist('hemisphere','var') || hemisphere==-1
    fwd_projection=@(lat, lon)ll2ps(lat, lon);
else
    fwd_projection=@(lat, lon)gl_ll2ps(lat, lon);
end
 
[~, j, k] = size(CP.wave_watts);

%---------------------calculate X (nadir) vector--------------------------%
[CP.xv_X, CP.yv_X, CP.zv_X]=ell2xyz(deg2rad(D.lat_20_ku(:)),deg2rad(D.lon_20_ku(:)),D.alt_20_ku(:),RefEll.a,RefEll.e2);
[temp.x, temp.y, temp.z]=ell2xyz(deg2rad(D.lat_20_ku(:)),deg2rad(D.lon_20_ku(:)),D.alt_20_ku(:)+1000,RefEll.a,RefEll.e2);
[CP.xv_X, CP.yv_X, CP.zv_X]=deal(CP.xv_X-temp.x, CP.yv_X-temp.y, CP.zv_X-temp.z);
Xnorm=sqrt(CP.xv_X.^2+ CP.yv_X.^2 + CP.zv_X.^2);
[CP.xv_X, CP.yv_X, CP.zv_X]=deal(CP.xv_X./Xnorm, CP.yv_X./Xnorm, CP.zv_X./Xnorm);

%--------------calculate V (spacecraft velocity) vector-------------------%
%Velocity vectors x,y,z in ITRF
CP.xv_V = D.sat_vel_vec_20_ku(1,:);
CP.yv_V = D.sat_vel_vec_20_ku(2,:);
CP.zv_V = D.sat_vel_vec_20_ku(3,:);

%normalize velocity vector
Vnorm = sqrt(CP.xv_V.^2 + CP.yv_V.^2 + CP.zv_V.^2);
CP.xv_V = CP.xv_V./Vnorm;
CP.yv_V = CP.yv_V./Vnorm;
CP.zv_V = CP.zv_V./Vnorm;


%-----------calculate Y (ellipsoidal velocity) and Z vector---------------%
X = [CP.xv_X(:) CP.yv_X(:) CP.zv_X(:)]'; %X vector
V = [CP.xv_V(:) CP.yv_V(:) CP.zv_V(:)]'; %V vector

%calculate Y vector (Cryosat-2 handbook eqn 2.1.3-1;
dotXV = dot(X,V);
dotXV = repmat(dotXV,[3 1]);
Y = V-(X.*dotXV);
% normalize:
Y=Y./repmat(sqrt(sum(Y.^2)), [3 1]);
CP.xv_Y = Y(1,:)';
CP.yv_Y = Y(2,:)';
CP.zv_Y = Y(3,:)';

%cross product of Y vector and X vector to find Z vector
CP.xv_Z = CP.yv_Y.*CP.zv_X - CP.zv_Y.*CP.yv_X;
CP.yv_Z = CP.zv_Y.*CP.xv_X - CP.xv_Y.*CP.zv_X;
CP.zv_Z = CP.xv_Y.*CP.yv_X - CP.yv_Y.*CP.xv_X;


%--------------------------determine ranges-------------------------------%


%sc position in CT system
[a,~,e2,~]=refell('WGS84'); %define CT system
D.lat_20_ku(D.lat_20_ku == 0) = NaN;
D.lon_20_ku(~isfinite(D.lat_20_ku)) = NaN;
[xp_sc,yp_sc,zp_sc]=ell2xyz(deg2rad(D.lat_20_ku(:)),deg2rad(D.lon_20_ku(:)), D.alt_20_ku(:),a,e2);
xp_sc = reshape(xp_sc,j,k);
yp_sc = reshape(yp_sc,j,k);
zp_sc = reshape(zp_sc,j,k);
[CP.xps_sc,CP.yps_sc]=fwd_projection(D.lat_20_ku(:),D.lon_20_ku(:)); %polar stereographic
CP.h_sc=reshape(D.alt_20_ku(:), j, k);
%calculate range
c = 299792458; %speed of light (m/s)
B = 3.2e8; %measured chirp bandwidth (Hz)
win_delay = D.window_del_20_ku.*(D.uso_cor_20_ku+1); %USO_Corr_factor, field 2;
Nsamps=size(D.power, 1);
range_surf_raw = (win_delay(burst_num)*c)/2 - (Nsamps*c)/(8*B) + ((sfc_bin-1)*c)/(4*B); %range(m)

sec_num=D.ind_meas_1hz_20_ku(burst_num)+1;
dry_trop = D.mod_dry_tropo_cor_01(sec_num);
wet_trop =  D.mod_wet_tropo_cor_01(sec_num);
model_ion =  D.iono_cor_01(sec_num);
ocean_loading_tide = D.load_tide_01(sec_num);
solidearth_tide = D.solid_earth_tide_01(sec_num);
geocentric_polar_tide =  D.pole_tide_01(sec_num);

range_surf = range_surf_raw(:)+dry_trop(:) +wet_trop(:) +model_ion(:) +ocean_loading_tide(:) +solidearth_tide(:) +geocentric_polar_tide(:);

% NOTE: baseline-D roll angle is in degrees
roll=D.off_nadir_roll_angle_str_20_ku*pi/180;

%-------------------find inferred angle from spacecraft-------------------%


lambda = 0.022084; %(m)
baseline = 1.1676; %(m)
 
[CP.lat_sc, CP.lon_sc]=deal(D.lat_20_ku(:), D.lon_20_ku(:));
CP.bad_file = 0;
delta_ph_ambig=2*pi*[-1 0 1];


% loop over three possible ambiguities
for ii = 1:3
    % 7/7/2016: The phase angles appear to need correction by a factor of .973 
    %(see Galin et al, 2013).
    %sfc_phase=D_in.phase*0.973;
    alpha = lambda*(sfc_phase+delta_ph_ambig(ii))/(2*pi*baseline); %alpha = inferred angle (rad) 512xjxk
    %alpha = repmat(alpha,[1 1 512]);
    %alpha = permute(alpha,[3 1 2]);
    R.inf_angle = (alpha+roll(burst_num)); %roll corrected inferred angle (rad)   !!! changed to + here!!!
    
    %calculate range vector
    R.xv_R = CP.xv_X(burst_num).*(range_surf.*cos(R.inf_angle))+CP.xv_Z(burst_num).*(range_surf.*sin(R.inf_angle));
    R.yv_R = CP.yv_X(burst_num).*(range_surf.*cos(R.inf_angle))+CP.yv_Z(burst_num).*(range_surf.*sin(R.inf_angle));
    R.zv_R = CP.zv_X(burst_num).*(range_surf.*cos(R.inf_angle))+CP.zv_Z(burst_num).*(range_surf.*sin(R.inf_angle));
    
    %normalize range vector
    Rnorm = sqrt(R.xv_R.^2 + R.yv_R.^2 + R.zv_R.^2);
    R.xv_R = R.xv_R./Rnorm;
    R.yv_R = R.yv_R./Rnorm;
    R.zv_R = R.zv_R./Rnorm;
    
    %---------------------calculate new point position------------------------%
    
    %new points in cartesian (m)
    xp_new = xp_sc(burst_num) + R.xv_R.*range_surf;
    yp_new = yp_sc(burst_num) + R.yv_R.*range_surf;
    zp_new = zp_sc(burst_num) + R.zv_R.*range_surf;
    
    %convert and resize new points to lat/lon/h (rad,rad,m)
    [R.lat_new, R.lon_new, R.h_new]=xyz2ell2(xp_new(:),yp_new(:),zp_new(:),a,e2);
    R.lat_new = rad2deg(R.lat_new);
    R.lon_new = rad2deg(R.lon_new);
    
    %convert and new and sc locations to polar sterographic. resize sc locations for plotting
    [R.xPS_new,R.yPS_new]=fwd_projection(R.lat_new,R.lon_new);    
    
    % export
    D_out.xPS(:,ii)=R.xPS_new(:)+1i*R.yPS_new;
    D_out.h(:,ii)=R.h_new;
    D_out.xv(:,ii)=R.xv_R(:)+1i*R.yv_R(:);
    D_out.zv(:,ii)=R.zv_R(:);
end

ff={'h_sc','xps_sc','yps_sc'};
for kf=1:length(ff)
    D_out.(ff{kf})=CP.(ff{kf})(burst_num);
end

D_out.range_surf=range_surf;

D_out.time=D.time_20_ku(burst_num)/24/3600+datenum('jan 1 2000');
D_out.time(D.time_20_ku==0)=NaN;
ff=fieldnames(D_in);
for kf=1:length(ff)
    D_out.(ff{kf})=D_in.(ff{kf});
end

XR=range(real(D_out.xPS(:)));
YR=range(imag(D_out.xPS(:)));
geoid=read_geotif_xy(params.geoid, XR+[-1e4, 1e4], YR+[-1e4, 1e4]);
if length(geoid.x)>1 && length(geoid.y)> 1
    D_out.geoid=interp2(geoid.x, geoid.y, geoid.z, real(D_out.xPS), imag(D_out.xPS));
else
    D_out.geoid=NaN(size(D_out.xPS));
end
    
    
bin_fn=round_to(D_out.xPS(isfinite(D_out.xPS)), 2.e4);
uBin=unique(bin_fn);
D_out.DEM=zeros(size(D_out.xPS))+NaN;
for kB=1:length(uBin)
    els=bin_fn==uBin(kB);
    XR=range(real(D_out.xPS(els)))+[-1e3, 1e3];
    YR=range(imag(D_out.xPS(els)))+[-1e3, 1e3];
    try
        DEM=read_geotif_xy(params.DEM, XR, YR);
        D_out.DEM(els)=interp2(DEM.x, DEM.y, DEM.z, real(D_out.xPS(els)), imag(D_out.xPS(els)));
    catch
        continue
    end
end
 