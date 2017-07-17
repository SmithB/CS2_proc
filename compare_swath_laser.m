function compare_swath_laser

out_file='laser_sw_comparison.mat';

ATM_master='/Volumes/insar7/gmap/oib_database/ATM_Qfit/Antarctica/master_index_h5.mat';
LVIS_master='/Volumes/insar6/gmap/oib_database/LVIS/Antarctica/master_index.mat';


load(ATM_master);
MI(1)=master_index;
load(LVIS_master);
MI(2)=master_index;
all_bins=cat(1, MI.bin_x)+1i*cat(1, MI.bin_y);
load SW_h5_C/south/AA_Shean_v1c_sigma4/master_index_h5.mat
I_sw=master_index;

bins=intersect(all_bins, I_sw.bin_x+1i*I_sw.bin_y);

[xg,yg]=meshgrid([-1 0 1]);
delta_bin=(xg(:)+1i*yg(:));

fields={'x','y','time','h', 'power','coherence','AD','error_composite', 'R_POCA', 'ambiguity','burst','abs_orbit', 'block_h_spread','count','phase','R_offnadir','dRange_POCA','POCA_flag'};

laser_fields={'x','y','z','time','N_50m','slope_x','slope_y'};

for kB=1:length(bins)
     
    D_sw=index_point_data_h5('read_from_index',  bins(kB)+ xg(:)+1i*yg(:),  I_sw, fields, []);
    
    D_sw=index_struct(D_sw, D_sw.coherence> 750  & D_sw.power > 1e-18 & D_sw.power < 1e-11 & D_sw.error_composite==0 & D_sw.count > 3 & D_sw.block_h_spread < 20);
    
  
    if isempty(D_sw.time); continue; end
    for kL=1:2;
        temp=index_point_data_h5('read_from_index', bins(kB)+delta_bin*1e4, MI(kL), laser_fields);
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
    D_sw.bin_fn=round_to(double(D_sw.x+1i*D_sw.y), 200);
    USB=unique(D_sw.bin_fn);
    USB=USB(ismember(USB, uL));
    USB1=unique(repmat(USB(:), [1, length(delta_bin(:))])+repmat(200*delta_bin(:).', [length(USB), 1]));
    DL=index_struct(DL, ismember(DL.bin_fn, USB1));
    D_sw=index_struct(D_sw, ismember(D_sw.bin_fn, USB1));
    if isempty(D_sw.time); continue; end
    
    sub_bins=unique(DL.bin_fn);
    for kS=1:length(sub_bins)
        els=ismember(DL.bin_fn, sub_bins(kS)+delta_bin*200);
        DL1=index_struct(DL, els);
        els_sw=D_sw.bin_fn==sub_bins(kS) & isfinite(D_sw.h);
        if ~any(els_sw); continue; end
        D_sw1=index_struct(D_sw, els_sw);
        
        D_sw1=index_struct(D_sw, els_sw);
        for kP=1:length(D_sw1.x);
            r=(DL1.x+1i*DL1.y)-(D_sw1.x(kP)+1i*D_sw1.y(kP));
            these=abs(r)< 100 & isfinite(DL1.z) & (abs(D_sw1.time(kP)-DL1.time) < 15);
            if ~any(these);
                continue
            end
            these=find(these);
            hmed=median(DL1.z(these));
            dh=abs(DL1.z(these)-hmed);
            this=these(find(dh==min(dh), 1, 'first'));
            D_sw1.h_laser_med(kP,1)=DL1.z(this);
            D_sw1.x_laser_med(kP,1)=DL1.x(this);
            D_sw1.y_laser_med(kP,1)=DL1.y(this);
            D_sw1.t_laser_med(kP,1)=DL1.time(this);
            D_sw1.slope_50m(kP,1)=DL1.slope_x(this)+1i*DL1.slope_y(this);
            D_sw1.N_50m(kP,1)=DL1.N_50m(this);
            D_sw1.h_laser_bar(kP,1)=mean(DL1.z(these));
            D_sw1.x_laser_bar(kP,1)=mean(DL1.x(these));
            D_sw1.y_laser_bar(kP,1)=mean(DL1.y(these));
            D_sw1.t_laser_bar(kP,1)=mean(DL1.time(these));
            D_sw1.h_laser_sigma(kP,1)=std(DL1.z(these));
            D_sw1.h_laser_iqr(kP,1)=iqr(DL1.z(these))/2;
            D_sw1.x_laser_sigma(kP,1)=sqrt(std(DL1.x(these)).^2+std(DL1.y(these)));
            D_sw1.N_laser(kP,1)=length(these);
        end
        if ~isfield(D_sw1,'N_laser'); continue; end
        good=find(D_sw1.N_laser >=1);
        if ~isempty(good)
            D_sw1=index_struct(D_sw1, good);
            if ~exist('D_both','var');
                D_both=D_sw1;
            else
                D_both(end+1)=D_sw1;
            end
        end
    end
end

f=fieldnames(D_both); for k=1:length(f); D_LCS.(f{k})=cat(1, D_both.(f{k})); end
D_LCS.year=(D_LCS.time-datenum('jan 1 2010'))/365.25+2010;
save(out_file,'D_both','D_LCS');

figure; 
uY=unique(round(D_LCS.year));
h=cheek_by_jowl(length(uY), 1, [0.1 0.1 0.8 0.8]); 
good=D_LCS.power > 1e-17 & D_LCS.power < 1e-13 & D_LCS.error_composite==0 & D_LCS.count > 3 & D_LCS.block_h_spread < 15;
for k=1:length(uY); 
    els=round(D_LCS.year)==uY(k) & D_LCS.h > 500 & good ;
    delta=D_LCS.h(els)-D_LCS.h_laser_med(els);
    axes(h(k)); 
    histogram(delta, -20:.05:20); 
    ht=text(10, mean(get(gca,'ylim')), {sprintf('med dh=%3.2f', median(delta)),  sprintf('iqr dh=%3.2f', iqr(delta)/2)});
end


