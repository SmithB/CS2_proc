function proc_L1_to_POCA_and_SW_BD(in_file, out_base, year)

try
    lat=read_nc_field(in_file,'/lat_20_ku');
catch
    fprintf(1,'proc_L1_to_POCA_and_SW: could not read latitude from file %s, returning\n', in_file);
    return
end

if isnumeric(year)
    year_str=num2str_year;
else
    year_str=year;
end


if max(lat(isfinite(lat))) < -40
    hemisphere=-1;
    out_dir=[out_base,'/AA/retrack/', year_str];
    params.DEM=[out_base,'/AA/DEM.tif'];
    params.geoid=[out_base,'/AA/geoid.tif'];
elseif min(lat(isfinite(lat)))>40
    hemisphere=1;
    out_dir=[out_base,'/GL/retrack/', year_str];
    params.DEM='GL/DEM.tif';
    params.geoid='GL/geoid.tif';
else
    disp('mid-latitude?  returning!')
    return
end

if ~exist(out_dir, 'dir')
    mkdir(out_dir)
end

[~, fname]=fileparts(in_file);

out_POCA=[out_dir,'/', fname,'.h5'];
out_SW=[out_dir,'/', fname,'_sw.h5'];
if exist(out_POCA,'file')
    return
end

[D_POCA, Dsw, orb]  = load_POCA_and_SW_BD(in_file, hemisphere, [], params);
if isempty(D_POCA); return; end
D_POCA.x=D_POCA.xPS; D_POCA=rmfield(D_POCA, 'xPS');
D_POCA.abs_orbit = repmat(orb.abs_orbit, size(D_POCA.x));
if ~isempty(Dsw)
    Dsw.x=Dsw.xPS;
    Dsw=rmfield(Dsw, 'xPS');
    Dsw.abs_orbit = repmat(orb.abs_orbit, size(Dsw.x));
    Dsw.AD=repmat(orb.AD, size(Dsw.x));
end

f = fieldnames(D_POCA);
for ii = 1:length(f)
    D_POCA.(f{ii}) = D_POCA.(f{ii})(:);
end

X0=unique(round_to(D_POCA.x, 1e4));
X0=X0(X0~=0 & isfinite(X0));

if ~isempty(D_POCA.h)
    [M, INDEX]=index_point_data_h5('build_index_and_sort', D_POCA, in_file, '/', 1e4, 10.1);
    if ~isempty(INDEX.bin_x)
        this_temp_file=[tempname,'.h5'];
        index_point_data_h5('write_h5', this_temp_file, M, INDEX, true);
        movefile(this_temp_file, out_POCA);
    end
end
if ~isempty(Dsw)
    X0=unique(round_to(Dsw.x, 1e4));
    X0=X0(X0~=0 & isfinite(X0));
    else
      fprintf(1, 'help! D_sw is empty for %s\n', out_SW)
end
fprintf(1, 'out_sw=%s\n', out_SW)
if ~isempty(Dsw) &&  ~isempty(Dsw.h);
    [M, INDEX]=index_point_data_h5('build_index_and_sort', Dsw, in_file, '/', 1e4, 10.3);
    if ~isempty(INDEX.bin_x)
        this_temp_file=[tempname,'.h5'];
        index_point_data_h5('write_h5', this_temp_file, M, INDEX, true);
        movefile(this_temp_file, out_SW);
    end
end



