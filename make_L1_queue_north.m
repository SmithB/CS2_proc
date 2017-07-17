function out_queue=make_L1_queue_north(in_dir, out_dir, DEM_file)

if ~exist('L1_DBL_list.txt','file')
    d=dir([in_dir,'/*.DBL']);
    fid=fopen('L1_DBL_list.txt','w');
    for k=1:length(d)
        fprintf(fid,'%s\n', d(k).name);
    end
    fclose(fid);
end
 
% get the list of files;     
file_list_FID=fopen('L1_DBL_list.txt','r');
 
if ~exist(out_dir,'dir'); mkdir(out_dir); end

out_queue={};
while ~feof(file_list_FID)
    in_file=deblank(fgetl(file_list_FID));
    
    [~, fname]=fileparts(in_file);
    out_POCA=[out_dir, fname,'.h5'];
    out_SW=[out_dir, fname,'_sw.h5'];
    if exist(out_POCA,'file')
        continue
    end
    out_queue{end+1}=sprintf('proc_L1_to_POCA_and_SW(''%s/%s'', ''%s'', 1, ''%s'');',in_dir, in_file, out_dir, DEM_file);
end

fclose(file_list_FID);

if nargout==0
    queue_list_name='L1_queue.txt';
    fid=fopen(queue_list_name,'w');  
    for k=1:length(out_queue);
        fprintf(fid, '%s\n', out_queue{k}); 
    end
    fclose(fid);
end