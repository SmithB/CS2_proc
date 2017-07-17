function compare_POCA_laser

ATM_master='/Volumes/insar7/gmap/oib_database/ATM_Qfit/Antarctica/master_index_h5.mat';
LVIS_master='/Volumes/insar6/gmap/oib_database/LVIS/Antarctica/master_index.mat';
out_file=sprintf('laser_POCA_comparison.mat');

load POCA_h5_C/south/AA_Shean_v1c_sigma4/master_index_h5.mat
I_POCA=master_index;

load(ATM_master);
MI(1)=master_index;
load(LVIS_master);
MI(2)=master_index; 

bins=intersect(cat(1, MI.bin_x)+1i*cat(1, MI.bin_y), I_POCA.bin_x+1i*I_POCA.bin_y);

[xg,yg]=meshgrid([-1 0 1]);
delta_bin=(xg(:)+1i*yg(:));
laser_fields={'x','y','z','time','N_50m','slope_x','slope_y'};
for kB=1:length(bins)
    %this_time_range=I_ATM.time_range{kB};
    D_POCA=index_point_data_h5('read_from_index', bins(kB), I_POCA,  {'x','y','h','time','AD','power','coherence','error_composite','abs_orbit'});
    % filter:
    D_POCA=index_struct(D_POCA, D_POCA.power> 3e-17 & D_POCA.coherence> 800 & D_POCA.error_composite==0 & D_POCA.power < 1e-12);
    %D_POCA=index_struct(D_POCA, ~ismember(D_POCA.abs_orbit, [ 1613 16900 16907 18117 18122, 18123, 24182 24183, 24185, 24190 24191 24192 24177 26305 26306]));
    %D_POCA=index_struct(D_POCA, D_POCA.time > min(this_time_range)-15 & D_POCA.time < max(this_time_range)+15);
    if isempty(D_POCA.time); continue; end
    for kL=1:2
        temp=index_point_data_h5('read_from_index', bins(kB)+delta_bin*1e4, MI(kL), {'x','y','z','time','N_50m','slope_x','slope_y'});
        if isfield(temp,'x')
            temp.bin_fn=round_to(double(temp.x+1i*temp.y), 200);
            DL(kL)=temp;
        else
            for kf=1:length(laser_fields); DL(kL).(laser_fields{kf})=[]; end; DL(kL).bin_fn=[];
        end
    end
    DL=cat_fields(DL(1), DL(2));
    
    uL=unique(DL.bin_fn(:));
    uL=unique(repmat(uL, [1, length(delta_bin(:))])+repmat(200*delta_bin(:).', [length(uL), 1]));
    
    uL=unique(repmat(uL, [1, length(delta_bin(:))])+repmat(200*delta_bin(:).', [length(uL), 1]));
    D_POCA.bin_fn=round_to(double(D_POCA.x+1i*D_POCA.y), 200);
    USB=unique(D_POCA.bin_fn);
    USB=USB(ismember(USB, uL));
    USB1=unique(repmat(USB(:), [1, length(delta_bin(:))])+repmat(200*delta_bin(:).', [length(USB), 1]));
    DL=index_struct(DL, ismember(DL.bin_fn, USB1));
    D_POCA=index_struct(D_POCA, ismember(D_POCA.bin_fn, USB1));
    if isempty(D_POCA.time); continue; end
    
    
    sub_bins=unique(DL.bin_fn);
    for kS=1:length(sub_bins)
        els=ismember(DL.bin_fn, sub_bins(kS)+delta_bin*200);
        DL1=index_struct(DL, els);
        els_sw=D_POCA.bin_fn==sub_bins(kS) & isfinite(D_POCA.h);
        if ~any(els_sw); continue; end
        D_POCA1=index_struct(D_POCA, els_sw);
        for kP=1:length(D_POCA1.x)
            r=(DL1.x+1i*DL1.y)-(D_POCA1.x(kP)+1i*D_POCA1.y(kP));
            these=abs(r)< 100 & isfinite(DL1.z) & (abs(D_POCA1.time(kP)-DL1.time) < 15);
            if ~any(these)
                continue
            end
            these=find(these);
            hmed=median(DL1.z(these));
            dh=abs(DL1.z(these)-hmed);
            this=these(find(dh==min(dh), 1, 'first'));
            D_POCA1.h_laser_med(kP,1)=DL1.z(this);
            D_POCA1.x_laser_med(kP,1)=DL1.x(this);
            D_POCA1.y_laser_med(kP,1)=DL1.y(this);
            D_POCA1.t_laser_med(kP,1)=DL1.time(this);
            D_POCA1.slope_50m(kP,1)=DL1.slope_x(this)+1i*DL1.slope_y(this);
            D_POCA1.N_50m(kP,1)=DL1.N_50m(this);
            D_POCA1.h_laser_bar(kP,1)=mean(DL1.z(these));
            D_POCA1.x_laser_bar(kP,1)=mean(DL1.x(these));
            D_POCA1.y_laser_bar(kP,1)=mean(DL1.y(these));
            D_POCA1.t_laser_bar(kP,1)=mean(DL1.time(these));
            D_POCA1.h_laser_sigma(kP,1)=std(DL1.z(these));
            D_POCA1.h_laser_iqr(kP,1)=iqr(DL1.z(these))/2;
            D_POCA1.x_laser_sigma(kP,1)=sqrt(std(DL1.x(these)).^2+std(DL1.y(these)));
            D_POCA1.N_laser(kP,1)=length(these);
        end
        if ~isfield(D_POCA1,'N_laser'); continue; end
        good=find(D_POCA1.N_laser >=1);
        if ~isempty(good)
            D_POCA1=index_struct(D_POCA1, good);
            if ~exist('D_both','var')
                D_both=D_POCA1;
            else
                D_both(end+1)=D_POCA1;
            end
        end
    end
end   
    
 
f=fieldnames(D_both); for k=1:length(f); D_LCS.(f{k})=cat(1, D_both.(f{k})); end
D_LCS.year=(D_LCS.time-datenum('jan 1 2010'))/365.25+2010;

save(out_file,'D_both','D_LCS');


figure; h=cheek_by_jowl(5, 1, [0.1 0.1 0.8 0.8]); 
 
good=D_LCS.power > 1e-16 &D_LCS .power < 1e-12 & D_LCS.error_composite==0;
uY=unique(round(D_LCS.year));
for k=1:5
    els=round(D_LCS.year)==uY(k) & D_LCS.h > 500 & good ;
    axes(h(k)); 
    histogram(D_LCS.h(els)-D_LCS.h_laser_med(els), -20:.05:20); 
    ht=text(10, mean(get(gca,'ylim')), sprintf('med dh=%3.2f, iqr dh=%3.2f', median(D_LCS.h(els)-D_LCS.h_laser_med(els)), iqr(D_LCS.h(els)-D_LCS.h_laser_med(els))/2)); 
end


