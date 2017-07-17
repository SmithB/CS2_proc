function proc_L1_to_POCA_and_SW(in_file, out_dir,  hemisphere, DEM)

if ~exist('hemisphere','var') || hemisphere==1
    DEM='/Volumes/insar2/gmap/gimp/gimp_90m_cubic.tif' ;
    hemisphere=1;
end

[~, fname]=fileparts(in_file);
%out_POCA=[out_dir,'/', fname,'.mat'];
%out_SW=[out_dir,'/', fname,'_sw.mat'];
out_POCA=[out_dir,'/', fname,'.h5'];
out_SW=[out_dir,'/', fname,'_sw.h5'];

if exist(out_POCA,'file')
    %return
end
%if lockfile_tool('lock', out_POCA)>0
%    return
%end

%fprintf(1, '%s\n', in_file);
%fprintf(1, '%s\n', out_POCA);
[D, Dsw, orb, out_table]  = read_CS2_file(in_file, hemisphere, [], struct('DEM', DEM));
if isempty(D); return; end
D.x=D.xPS; D=rmfield(D, 'xPS');
D.abs_orbit = repmat(orb.abs_orbit, size(D.x));
if ~isempty(Dsw)
    Dsw.x=Dsw.xPS;
    Dsw=rmfield(Dsw, 'xPS');
    [~, ind]=ismember(Dsw.burst, out_table.burst);
    copy_fields={'AD','error_composite'};
    for kf=1:length(copy_fields)
        Dsw.(copy_fields{kf})=out_table.(copy_fields{kf})(ind);
    end
    Dsw.abs_orbit = repmat(orb.abs_orbit, size(Dsw.x));
end

f = fieldnames(D);
for ii = 1:length(f)
    D.(f{ii}) = D.(f{ii})(:);
end

X0=unique(round_to(D.x, 1e4));
X0=X0(X0~=0 & isfinite(X0));


%save(out_POCA, 'D','X0');
if ~isempty(D.h);
    [M, INDEX]=index_point_data_h5('build_index_and_sort', D, in_file, '/', 1e4, 10.1);
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
%save( out_SW, 'Dsw','X0')

%lockfile_tool('unlock', out_POCA);




