function [out_queue, status, INDEX]=check_L1_queue_north(in_dir, out_dir, DEM_file)

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
[~, out]=unix('wc -l L1_DBL_list.txt | awk ''{print $1}''');
N_in=str2num(out);
INDEX=repmat(struct('POCA', struct('bin_x', [],'bin_y', []),'SW', struct('bin_x', [],'bin_y', [])), [N_in, 1]);
status=repmat(struct('POCA', 0,'SW', 0), [N_in, 1]); 

for count=1:N_in
    in_file=deblank(fgetl(file_list_FID));
    [~, fname]=fileparts(in_file);
    out_POCA=[out_dir,'/' fname,'.h5'];
    out_SW=[out_dir, '/', fname,'_sw.h5'];
    if ~exist(out_POCA,'file')
        status(count).POCA=1;
    else
        try
            INDEX(count).POCA.bin_x=h5read(out_POCA, '/INDEX/bin_x');
            INDEX(count).POCA.bin_y=h5read(out_POCA, '/INDEX/bin_y');
            INDEX(count).POCA.file_num=count;
        catch
            status(count).POCA=2;
        end
    end
    if ~exist(out_SW,'file');
        status(count).SW=1;
    else
        try
            INDEX(count).SW.bin_x=h5read(out_SW, '/INDEX/bin_x');
            INDEX(count).SW.bin_y=h5read(out_SW, '/INDEX/bin_y');
        catch
            status(count).SW=2;
        end
    end
    if status(count).POCA>0 && status(count).SW>0;
        out_queue{count}=sprintf('proc_L1_to_POCA_and_SW(''%s/%s'', ''%s'', 1, ''%s'');',in_dir, in_file, out_dir, DEM_file);
    end
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



% if there are problems with some files that dont' have good 'x' fields:


thedir='/Volumes/insar6/ben/Cryosat/slope_picks_C/north/gimp_v1c_sigma4/';
[~, files]=unix(sprintf('ls %s | grep h5', thedir));
files=strsplit(deblank(files));

bad_files={};
for k=1:length(files)
    try
        ii=h5info([thedir,'/',files{k}],'/x');
    catch
        disp(files{k});
        bad_files{end+1}=files{k};
    end
end

% cleanup:
if 0
    for k=1:length(bad_files);
        temp=strrep(strrep(bad_files{k},'_sw.h5',''),'.h5',''); 
        [ss, theline]=unix(sprintf('grep %s L1_queue.txt', temp));
        if ss==0;
            unix(sprintf('rm %s/%s*', thedir, temp));
            eval(theline); 
        end
    end
end


