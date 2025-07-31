[~, host]=unix('hostname');
host=deblank(host);

setappdata(0,'hemisphere', -1);
info=default_info;
addpath([info.data_dir,'/../Tides/']);


thename='AA_Helm_v2c_sigma4';
%DEM_file='/Volumes/ice1/ben/scale_maps/wais_coastline_mos_20151013_full_embedded_bedmap2_width3px-tile-0_100m.tif';
DEM_file='/Volumes/ice1/ben/helm_DEM/AA/helm_DEM_with_EGM2008.tif';%
    
L1_dir='/Volumes/insar6/ben/Cryosat/L1_C/south/';
L2_dir='/Volumes/insar6/ben/Cryosat/';
slope_picks_dir=sprintf('%s/slope_picks_C/south/%s/', L2_dir, thename);
h5_dir_sw=sprintf('%s/SW_h5_C/south/%s/', L2_dir, thename);
h5_dir_poca=sprintf('%s/POCA_h5_C/south/%s/', L2_dir, thename);

if false
    temp=gcp;
    if temp.NumWorkers ~=12
        delete(gcp)
        parpool(12)
    end
    % run make_L1_queue_south
    L=make_L1_queue_south(L1_dir, slope_picks_dir, DEM_file);
    % now run it using parallel toolbox
    
    
    if strcmp(host,'store')
        FirstJob=ceil(length(L)/2+.1); LastJob=length(L);
    else
        FirstJob=1; LastJob=floor(length(L)/2);
    end
    
    %parfor k=FirstJob:LastJob
    parfor k=1:length(L)
        pareval(L{k});
    end
end


if true
    if false
        if exist(h5_dir_poca,'dir')
            delete([h5_dir_sw,'/*']);
            delete([h5_dir_poca,'/*']);
        else
            mkdir(h5_dir_poca);
            mkdir(h5_dir_sw);
        end
    end
    queue_sw=make_AA_h5('make_queue', slope_picks_dir, h5_dir_sw, 5e4, 'sw');
    queue_POCA=make_AA_h5('make_queue', slope_picks_dir, h5_dir_poca, 2e5, 'POCA');
    
    if false
        disp('now running the POCA queues in parallel')
        delete(gcp('nocreate'));
        
        if strcmp(host,'store')
            parpool(8);
            FirstJob=ceil(length(queue_POCA)/2+.1); LastJob=length(queue_POCA);
        else
            parpool(4);
            FirstJob=1; LastJob=floor(length(queue_POCA)/2);
        end
        parfor k=FirstJob:LastJob
            pareval(queue_POCA{k});
        end
    end
    make_AA_h5('make_index',  h5_dir_poca);
    %compare_POCA_laser
    
    disp('now running the swath queues in parallel')
    delete(gcp('nocreate'));
    
    if strcmp(host,'store')
        parpool(2);
        FirstJob=ceil(length(queue_sw)/2+.1); LastJob=length(queue_sw);
    else
        parpool(2);
        FirstJob=1; LastJob=floor(length(queue_sw)/2);
    end
    fprintf(1, '%d jobs in queue', length(queue_sw));
    for k=1:length(queue_sw);
    %parfor k=FirstJob:LastJob
        eval(queue_sw{k});
    end
    make_AA_h5('make_index',  h5_dir_sw);
    % compare_swath_laser
end