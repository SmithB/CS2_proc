

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