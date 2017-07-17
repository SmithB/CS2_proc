function  varargout=make_CS_h5(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end


function queue_list=make_queue(  master_index_file, out_dir, dx0, mode)

load(master_index_file,'master_index');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
uBins=unique(round_to( master_index.bin_x+1i*master_index.bin_y, dx0));
fid=fopen([out_dir,'/bin_list.txt'],'w');
for k=1:length(uBins)
    fprintf(fid,'%s_%dE_%dN.h5\n', mode, real(uBins(k))/1000, imag(uBins(k))/1000);
end
fclose(fid);
 
fid_list=fopen([out_dir,'/bin_list.txt'],'r');

count=0;
while ~feof(fid_list)
    fname=deblank(fgetl(fid_list));
    count=count+1;
    queue_list{count}=sprintf( 'make_CS_h5(''run_bin'', ''%s'', ''%s'', ''%s'', %d);', master_index_file, out_dir, fname, dx0);
end
fclose(fid_list);

if nargout==0
    fid_queue=fopen(sprintf('queue_for_CS%s.txt',mode),'w');
    for k=1:length(queue_list)
        fprintf(fid_queue, '%s\n', queue_list{k});
    end
    fclose(fid_queue);
end

%-----------------------------------------------------------------
function run_bin(master_index_file, out_dir, out_file, dx0)

master_index=getappdata(0,'master_index'); 
if isempty(master_index)
    load(master_index_file);
    setappdata(0,'master_index', master_index);
end

fields=strsplit(h5readatt([master_index.base,'/',master_index.out_file{1}], '/','fieldname_list'));

fprintf(1,'make_CS_h5: working on %s\n', out_file);

temp=regexp(out_file,'([-+]*\d+)E_([-+]*\d+)N','tokens');
x0=str2double(temp{1}{1})*1000;
y0=str2double(temp{1}{2})*1000;

all_bins=master_index.bin_x+1i*master_index.bin_y;
these_bins=all_bins(round_to(all_bins, dx0)==(x0+1i*y0));
if isempty(these_bins)
    return
end

D=index_point_data_h5('read_from_index', these_bins, master_index, fields);
 
these=D.x>x0-dx0/2 & D.x <= x0+dx0/2  & D.y>y0-dx0/2 & D.y <= y0+dx0/2 & isfinite(D.x);
if ~any(these); return; end
if ~all(these);
D=index_struct(D, these);
end

% check this
[M, INDEX]=index_point_data_h5('build_index_and_sort', D, [], out_dir, 1e4, 10.2);

index_point_data_h5('write_h5',[out_dir,'/',out_file], M, INDEX, true);
%save([out_dir,'/', strrep(out_file,'.h5','_index.mat')], 'INDEX');

%-----------------------------------------------------------------
function make_index(out_dir)
 
% make the master index
[~, h5_files]=unix(['ls ', out_dir,'/*.h5']);
h5_files=strsplit(deblank(h5_files));
for k=1:length(h5_files);
    try
        II=h5info(h5_files{k},'/INDEX');
        disp(II.Filename);
    catch
        fprintf(1,'failure for %s\n', h5_files{k});
        delete(h5_files{k})
    end
end

INDEX_list=index_point_data_h5('collect_index', {out_dir}, '.h5',  {'./'});

master_index=index_point_data_h5('make_master_index', INDEX_list, [out_dir,'/master_index_h5.mat'],  out_dir);
