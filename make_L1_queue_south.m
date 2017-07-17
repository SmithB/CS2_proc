
function out_queue=make_L1_queue_south(in_dir, out_dir, DEM_file, file_list)

 
% get the list of files;     
file_list_FID=fopen('../file_list.txt','r');
 
if ~exist(out_dir,'dir'); mkdir(out_dir); end

if ~exist('file_list','var')
    [~, out]=unix(sprintf('ls %s | grep DBL', in_dir));
    file_list=strsplit(deblank(out));
end

out_queue={};
for in_file=file_list    
    [~, fname]=fileparts(in_file{1});
    out_POCA=[out_dir, fname,'.mat'];
    out_SW=[out_dir, fname,'_sw.mat'];
    if exist(out_POCA,'file')
        continue
    end
    out_queue{end+1}=sprintf('proc_L1_to_POCA_and_SW(''%s/%s'', ''%s'', -1, ''%s'');',in_dir, in_file{1}, out_dir, DEM_file);
end

 
if nargout==0
    queue_list_name='L1_queue.txt';
    fid=fopen(queue_list_name,'w');  
    for k=1:length(out_queue);
        fprintf(fid, '%s\n', out{k}); 
    end
    fclose(fid);
end