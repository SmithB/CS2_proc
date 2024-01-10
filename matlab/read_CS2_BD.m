function [D, orb]=read_CS2_BD(file)

fields={'alt_20_ku', 'coherence_waveform_20_ku', 'iono_cor_01', ...
    'ind_meas_1hz_20_ku', ...
    'lat_20_ku', ...
    'lon_20_ku', 'load_tide_01',...
    'flag_mcd_20_ku', 'mod_dry_tropo_cor_01', 'mod_wet_tropo_cor_01', ...
    'noise_power_20_ku', 'ocean_tide_eq_01', 'off_nadir_roll_angle_str_20_ku',...
    'ph_diff_waveform_20_ku', 'pole_tide_01', 'sat_vel_vec_20_ku',...
    'rec_count_20_ku', 'solid_earth_tide_01', 'time_20_ku', 'transmit_pwr_20_ku', 'uso_cor_20_ku' ...
    'window_del_20_ku'};
%?iono_cor_gim?    

D=struct();
for k=1:length(fields)
    D.(fields{k})=read_nc_field(file, fields{k});
end

orb.AD=((h5readatt(file,'/','ascending_flag')=='A')-0.5)*2;
orb.abs_orbit=ncreadatt(file,'/', 'abs_orbit_start');
orb.rel_orbit=ncreadatt(file,'/', 'rel_orbit_number');
orb.cycle=ncreadatt(file,'/', 'cycle_number');


% from the netcdf file description for echo_scale_factor_20_ku
% The 20Hz power waveform scaling factor, computed in order to best fit each 
% waveform within 2 bytes. The scaling, needed to convert the L1B 1Hz average 
% power waveform into Watts, is applied as follows: 
%  pwr_waveform_20_ku(time_20_ku,ns_20_ku)*echo_scale_factor_20_ku(time_20_ku)*2^echo_scale_pwr_20_ku(time_20_ku).";
N_samps=size(D.coherence_waveform_20_ku, 1);
D.power=read_nc_field(file, 'pwr_waveform_20_ku') .* ...
    repmat(read_nc_field(file, 'echo_scale_factor_20_ku')', [N_samps,1]) .*...
    2.^repmat(read_nc_field(file, 'echo_scale_pwr_20_ku')', [N_samps,1]);




