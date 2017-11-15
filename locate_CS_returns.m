function D_out=locate_CS_returns(D_in, L1b, CP, hemisphere)

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
[CP.xv_X, CP.yv_X, CP.zv_X]=ell2xyz(deg2rad(L1b.GEO.LAT(:)),deg2rad(L1b.GEO.LON(:)),L1b.GEO.H(:),RefEll.a,RefEll.e2);
[temp.x, temp.y, temp.z]=ell2xyz(deg2rad(L1b.GEO.LAT(:)),deg2rad(L1b.GEO.LON(:)),L1b.GEO.H(:)+1000,RefEll.a,RefEll.e2);
[CP.xv_X, CP.yv_X, CP.zv_X]=deal(CP.xv_X-temp.x, CP.yv_X-temp.y, CP.zv_X-temp.z);
Xnorm=sqrt(CP.xv_X.^2+ CP.yv_X.^2 + CP.zv_X.^2);
[CP.xv_X, CP.yv_X, CP.zv_X]=deal(CP.xv_X./Xnorm, CP.yv_X./Xnorm, CP.zv_X./Xnorm);

%--------------calculate V (spacecraft velocity) vector-------------------%
%Velocity vectors x,y,z in ITRF
CP.xv_V = L1b.GEO.V.Vx(:);
CP.yv_V = L1b.GEO.V.Vy(:);
CP.zv_V = L1b.GEO.V.Vz(:);

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

% %normalize Y vector
% bvnorm = sqrt(CP.xv_Y.^2 + CP.yv_Y.^2 + CP.zv_Y.^2); bvnorm=bvnorm(:);
% CP.xv_Y = CP.xv_Y./bvnorm;
% CP.yv_Y = CP.yv_Y./bvnorm;
% CP.zv_Y = CP.zv_Y./bvnorm;
 
%cross product of Y vector and X vector to find Z vector
CP.xv_Z = CP.yv_Y.*CP.zv_X - CP.zv_Y.*CP.yv_X;
CP.yv_Z = CP.zv_Y.*CP.xv_X - CP.xv_Y.*CP.zv_X;
CP.zv_Z = CP.xv_Y.*CP.yv_X - CP.yv_Y.*CP.xv_X;


%--------------------------determine ranges-------------------------------%


%sc position in CT system
[a,b,e2,finv]=refell('WGS84'); %define CT system
L1b.GEO.LAT(L1b.GEO.LAT == 0) = NaN;
L1b.GEO.LON(L1b.GEO.LON == 0) = NaN;
[xp_sc,yp_sc,zp_sc]=ell2xyz(deg2rad(L1b.GEO.LAT(:)),deg2rad(L1b.GEO.LON(:)),L1b.GEO.H(:),a,e2);
xp_sc = reshape(xp_sc,j,k);
yp_sc = reshape(yp_sc,j,k);
CP.h_sc = reshape(zp_sc,j,k);
[CP.xps_sc,CP.yps_sc]=fwd_projection(L1b.GEO.LAT(:),L1b.GEO.LON(:)); %polar stereographic

%calculate range
c = 299792458; %speed of light (m/s)
B = 3.2e8; %measured chirp bandwidth (Hz)
win_delay = L1b.MEA.win_delay.*(L1b.GEO.USO+1); %USO_Corr_factor, field 2;
range_surf = (win_delay(burst_num)*c)/2 - (size(L1b.SIN.data,1)*c)/(8*B) + ((sfc_bin-1)*c)/(4*B); %range(m)

%atmospheric and tidal corrections to range
% OLD VERSION (as of 5/15/2016)
% dry_trop = repmat(L1b.COR.dry_trop,[j 1]);
% wet_trop = repmat(L1b.COR.wet_trop,[j 1]);
% model_ion = repmat(L1b.COR.model_ion,[j 1]);
% ocean_loading_tide = repmat(L1b.COR.ocean_loading_tide,[j 1]);
% solidearth_tide = repmat(L1b.COR.solidearth_tide,[j 1]);
% geocentric_polar_tide = repmat(L1b.COR.geocentric_polar_tide,[j 1]);
%range_surf = range_surf+dry_trop(burst_num)+wet_trop(burst_num)+model_ion(burst_num)+ocean_loading_tide(burst_num)+solidearth_tide(burst_num)+geocentric_polar_tide(burst_num);


% NEW VERSION (as of 5/15/2016)
sec_num=CP.sec_number(burst_num);
dry_trop = L1b.COR.dry_trop(sec_num);
wet_trop =  L1b.COR.wet_trop(sec_num);
model_ion =  L1b.COR.model_ion(sec_num);
ocean_loading_tide = L1b.COR.ocean_loading_tide(sec_num);
solidearth_tide = L1b.COR.solidearth_tide(sec_num);
geocentric_polar_tide =  L1b.COR.geocentric_polar_tide(sec_num);


range_surf = range_surf(:)+dry_trop(:) +wet_trop(:) +model_ion(:) +ocean_loading_tide(:) +solidearth_tide(:) +geocentric_polar_tide(:);


%-------------------find inferred angle from spacecraft-------------------%

if isfield(L1b.GEO, 'Antenna_Bench_Roll');
    roll=L1b.GEO.Antenna_Bench_Roll*pi/180;
else
    roll=+L1b.GEO.BaseLine.X; %roll correction (rad)
end

% Attempt to read the roll correction:
burst_time=L1b.GEO.Serial_Sec_Num/24/3600+datenum('jan 1 2000');
roll_corrected=read_roll_correction(burst_time)*pi/180;
delta_roll=roll_corrected-roll;

% if the roll correction isn't available, subtract 0.0075 deg
bad_RC=~isfinite(roll_corrected);
roll_corrected(bad_RC)=roll(bad_RC)-.0075*pi/180;

roll(isfinite(roll_corrected))=roll_corrected(isfinite(roll_corrected));

lambda = 0.022084; %(m)
baseline = 1.1676; %(m)
 
[CP.lat_sc, CP.lon_sc]=deal(L1b.GEO.LAT(:), L1b.GEO.LON(:));
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
    zp_new = CP.h_sc(burst_num) + R.zv_R.*range_surf;
    
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

D_out.time=L1b.GEO.Serial_Sec_Num(burst_num)/24/3600+datenum('jan 1 2000');
D_out.time(L1b.GEO.Serial_Sec_Num(burst_num)==0)=NaN;
ff=fieldnames(D_in);
for kf=1:length(ff)
    D_out.(ff{kf})=D_in.(ff{kf});
end

D_out.delta_roll=delta_roll(burst_num);
 