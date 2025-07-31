[~, host]=unix('hostname');
host=deblank(host);

setappdata(0,'hemisphere', 1);
info=default_info;
addpath([info.data_dir,'/../Tides/']);

thename='gimp_v1c_sigma4';
DEM_file='/Volumes/insar9/ben/GIMP2_DEM/v2fit_120.tif';

L1_dir='/Volumes/insar6/ben/Cryosat/L1_C/north';
L2_dir='/Volumes/insar6/ben/Cryosat/';
slope_picks_dir=sprintf('%s/slope_picks_C/north/%s', L2_dir, thename);
h5_dir_sw=sprintf('%s/SW_h5_C/north/%s', L2_dir, thename);
h5_dir_poca=sprintf('%s/POCA_h5_C/north/%s', L2_dir, thename);

 
if true
    ! rm -rf par_runs
    make_L1_queue_north(L1_dir, slope_picks_dir, DEM_file);
    % start a new tmux session (if needed)
    ! bash -c "tmux info || (tmux new-session -d)"
    % start the parallel boss
    ! tmux new-window "parallel_boss.py -m L1_queue.txt ; parallel_boss.py"
    ! run_workers_in_tmux -m 5
    % Now run around to the other machines and start a bunch of other workers.  OR don't, if this is an easy job.
    pause(30);
    while unix('ls par_run/running')==0
        pause(60);
    end
    ! rm -rf par_runs
end


if true
    if true
        if exist(h5_dir_poca,'dir')
            delete([h5_dir_sw,'/*']);
            delete([h5_dir_poca,'/*']);
        else
            mkdir(h5_dir_poca);
            mkdir(h5_dir_sw);
        end
    end
    
    INDEX_list_POCA=index_point_data_h5('collect_index', {slope_picks_dir}, '.*C....h5',{'./'});
    % there's a problem with the time ranges in the indices-- fix:
    for k=1:length(INDEX_list_POCA)
        if size(INDEX_list_POCA(k).time_range, 2)==1 && size(INDEX_list_POCA(k).time_range, 1)==2*size(INDEX_list_POCA(k).bin_x,1) 
            INDEX_list_POCA(k).time_range=reshape(INDEX_list_POCA(k).time_range, [size(INDEX_list_POCA(k).bin_x,1), 2]); 
        end
    end
    [master_index_POCA, INDEX_POCA]=index_point_data_h5('make_master_index', INDEX_list_POCA, [slope_picks_dir,'/master_index_POCA.mat'], slope_picks_dir);
    
    INDEX_list_SW=index_point_data_h5('collect_index', {slope_picks_dir}, '.*C..._sw.h5',{'./'});
    % there's a problem with the time ranges in the indices-- fix:
    for k=1:length(INDEX_list_SW)
        if size(INDEX_list_SW(k).time_range, 2)==1 && size(INDEX_list_SW(k).time_range, 1)==2*size(INDEX_list_SW(k).bin_x,1) 
            INDEX_list_SW(k).time_range=reshape(INDEX_list_SW(k).time_range, [size(INDEX_list_SW(k).bin_x,1), 2]); 
        end
    end
    [master_index_SW, INDEX_SW]=index_point_data_h5('make_master_index', INDEX_list_SW, [slope_picks_dir,'/master_index_SW.mat'], slope_picks_dir);

    
    make_CS_h5('make_queue',  [slope_picks_dir,'/master_index_SW.mat'], h5_dir_sw,  5e4, 'sw');
    make_CS_h5('make_index',  h5_dir_sw);    
 
    make_CS_h5('make_queue', [slope_picks_dir,'/master_index_POCA.mat'], h5_dir_poca, 2e5, 'POCA');
    make_CS_h5('make_index',  h5_dir_poca);
    
     
    % compare_swath_laser
end