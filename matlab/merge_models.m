function [ dZ, z0, bias, Flist, x0, seasonal_cycle, bias_list]=merge_models(dirname, XR, YR);
% [ dZ, z0, bias, Flist, x0, seasonal_cycle, bias_list]=merge_models(dirname, XR, YR);
% merges Cryosat-2 models in directory 'dirname', optionally producing
% limiting the output range to XR and YR

files=dir([dirname,'/fit*.mat']);

for k=1:length(files)
    temp=regexp(files(k).name,'fit_(.*)_(.*).mat','tokens'); XY(k,:)=str2num([temp{1}{1},' ', temp{1}{2}])*1000; 
end
dX=median([diff(unique(XY(:,1))); diff(unique(XY(:,1)))]) ;
if exist('XR','var');
    these_files=XY(:,1) >= XR(1)-dX & XY(:,1) <= XR(2)+dX & XY(:,2) >= YR(1)-dX & XY(:,2) <= YR(2)+dX;
    files=files(these_files);
end
good=false(length(files),1);
for kf=1:length(files)
    
    try
        temp=load([dirname,'/',files(kf).name],'M');
        if ~(any(temp.M.YR>=YR(1) & temp.M.YR <=YR(2)) &&...
                any(temp.M.XR>=XR(1) &  temp.M.XR<=XR(2)))
            continue
        end
        L(kf)=load([dirname,'/',files(kf).name],'S','M');
        
    catch
        continue
    end
    if isstruct(L(kf).S)   
        temp=load([dirname,'/',files(kf).name],'D');
        [CC(kf).swath, CC(kf).poca]=count_swath_poca(index_struct(temp.D, temp.D.good), L(kf).S.grids);           
        XR_sub(kf,:)=range(L(kf).S.grids.dz.ctrs{2});
        YR_sub(kf,:)=range(L(kf).S.grids.dz.ctrs{1});
        x0(kf)=mean(L(kf).S.grids.dz.ctrs{2})+1i*mean(L(kf).S.grids.dz.ctrs{1});
        seasonal_cycle(kf,:)=L(kf).S.seasonal_cycle;
        xw(kf)=diff(range(L(kf).S.grids.dz.ctrs{2}));
        Flist{kf}=[dirname,'/', files(kf).name];
        good(kf)=true;
    end
end
XR_sub=XR_sub(good,:); 
YR_sub=YR_sub(good,:); 
seasonal_cycle=seasonal_cycle(good,:);
xw=xw(good); 
Flist=Flist(good);
L=L(good);
x0=x0(good);

dx_ctr=median(diff(unique(real(x0))));
overlap=(median(xw)-dx_ctr);
N_clip=floor((overlap/L(1).M.dx)/2);
dx=L(1).M.dx;

%if ~exist('XR','var');
XR=[min(XR_sub(:,1))+ceil(N_clip/2)*dx max(XR_sub(:,2))-ceil(N_clip/2)*dx];
YR=[min(YR_sub(:,1))+ceil(N_clip/2)*dx max(YR_sub(:,2))-ceil(N_clip/2)*dx];
%end
HW=hanning(N_clip).^2; HW=HW/sum(HW);

dZ.x=XR(1):dx:XR(2);
dZ.y=YR(1):dx:YR(2);
dZ.t=L(1).S.grids.dz.ctrs{3};
dZ.z=NaN(length(dZ.y), length(dZ.x), length(dZ.t));
dZ.count=dZ.z(:,:,1);
dZ.bias=dZ.count;
dZ.season_count=dZ.z;
bias.z=dZ.z(:,:,1);
dZ.sigma=zeros(size(dZ.z));


[zw, sw, scw, w]=deal(zeros(size(dZ.z)));
[ bw, cw]=deal(zeros(size(dZ.z(:,:,1))));
for kf=1:length(L);
    if ~isstruct(L(kf).S);
        continue
    end
    sz=size(L(kf).S.dz);
    sub_rows=ceil(N_clip+1):sz(1)-floor(N_clip+1);
    sub_cols=ceil(N_clip+1):sz(2)-floor(N_clip+1);   
    
    [~, sub_row_ind, ~]=intersect(L(kf).S.grids.dz.ctrs{1}(sub_rows), dZ.y);
    [~, sub_col_ind, ~]=intersect(L(kf).S.grids.dz.ctrs{2}(sub_cols), dZ.x);
    
    temp=zeros(size(L(kf).S.dz(:,:,1)));
    temp(sub_rows(sub_row_ind), sub_cols(sub_col_ind))=1;
    temp=conv2(conv2(temp, HW,'same'), HW','same');
    
    [~, sub_row_ind, rows0]=intersect(L(kf).S.grids.dz.ctrs{1}, dZ.y);
    [~, sub_col_ind, cols0]=intersect(L(kf).S.grids.dz.ctrs{2}, dZ.x);
    
    zw(rows0, cols0,:)=zw(rows0, cols0,:)+repmat(temp(sub_row_ind, sub_col_ind), [1 1 size(L(kf).S.dz,3)]).*L(kf).S.dz(sub_row_ind, sub_col_ind,:);
    w(rows0, cols0,:)=w(rows0, cols0,:)+repmat(temp(sub_row_ind, sub_col_ind), [1 1 size(L(kf).S.dz,3)]);
    scw(rows0, cols0,:)=scw(rows0, cols0,:)+repmat(temp(sub_row_ind, sub_col_ind), [1 1 size(L(kf).S.dz,3)]).*L(kf).S.count_by_season(sub_row_ind, sub_col_ind,:);
   
    %dZ.z(rows0, cols0, :)=L(kf).S.dz(sub_rows(sub_row_ind), sub_cols(sub_col_ind), :);    
    %bias.z(rows0, cols0)=L(kf).S.bias(sub_rows(sub_row_ind), sub_cols(sub_col_ind));
    bw(rows0, cols0)=bw(rows0, cols0)+ temp(sub_row_ind, sub_col_ind).*L(kf).S.bias(sub_row_ind, sub_col_ind);
    
    cw(rows0, cols0)=cw(rows0, cols0)+temp(sub_row_ind, sub_col_ind).*L(kf).S.count(sub_row_ind, sub_col_ind); 
    if isfield(L(kf).S,'sigma_dz')
        sw(rows0, cols0, :)=sw(rows0, cols0,:)+repmat(temp(sub_row_ind, sub_col_ind), [1 1 size(L(kf).S.dz,3)]).*L(kf).S.sigma_dz(sub_row_ind, sub_col_ind,:);
    end
    
    if isfield(L(kf).S,'roll_bias');
        bias_list.roll(kf)=L(kf).S.roll_bias;
    end
    if isfield(L(kf).S,'perswath_bias');
        bias_list.perswath(kf)=L(kf).S.perswath_bias;
    end
    
end

els1=w(:,:,1)>0.001;
els=w>0.001;
dZ.z(els)=zw(els)./w(els);
dZ.sigma(els)=sw(els)./w(els);
w1=w(:,:,1);
dZ.count(els1)=cw(els1)./w1(els1);
dZ.bias(els1)=bw(els1)./w1(els1);
dZ.season_count(els)=scw(els)./w(els);

dx=L(1).M.dx0;
N_clip=floor((overlap/dx)/2);
HW=hanning(N_clip).^2; HW=HW/sum(HW);

z0.x=XR(1):dx:XR(2);
z0.y=YR(1):dx:YR(2);
[z0.z, z0.count_POCA, z0.count_swath]=deal(NaN(length(z0.y), length(z0.x)));
[zw, w, cpw, csw]=deal(zeros(size(z0.z)));

for kf=1:length(L);
    sz=size(L(kf).S.z0);
    sub_rows=ceil(N_clip+1):sz(1)-floor(N_clip+1);
    sub_cols=ceil(N_clip+1):sz(2)-floor(N_clip+1);
    [~, sub_row_ind]=intersect(L(kf).S.grids.z0.ctrs{1}(sub_cols), z0.y);
    [~, sub_col_ind]=intersect(L(kf).S.grids.z0.ctrs{2}(sub_rows), z0.x);
    
    [~, S_row_ind, r0]=intersect(L(kf).S.grids.z0.ctrs{1}, z0.y);
    [~, S_col_ind, c0]=intersect(L(kf).S.grids.z0.ctrs{2}, z0.x);
    
    temp=zeros(size(L(kf).S.z0));
    temp(sub_rows(sub_row_ind), sub_cols(sub_col_ind))=1;
    temp=conv2(conv2(temp, HW,'same'), HW','same');

    w(r0, c0)=w(r0, c0)+temp(S_row_ind, S_col_ind);
    zw(r0, c0)=zw(r0, c0)+temp(S_row_ind, S_col_ind).*L(kf).S.z0(S_row_ind, S_col_ind);
    if ~isempty(CC(kf).poca)
        cpw(r0, c0)=cpw(r0, c0)+temp(S_row_ind, S_col_ind).*CC(kf).poca(S_row_ind, S_col_ind);
        csw(r0, c0)=csw(r0, c0)+temp(S_row_ind, S_col_ind).*CC(kf).swath(S_row_ind, S_col_ind);
    end
end
z0.z(w>0.001)=zw(w>0.001)./w(w>0.001);
z0.count_POCA(w>0.001)=cpw(w>0.001)./w(w>0.001);
z0.count_swath(w>0.001)=csw(w>0.001)./w(w>0.001);

if ~exist('bias_list','var'); bias_list=[]; end

if 0
h=cheek_by_jowl(5,2, [0 0 1 1])';
for k=1:size(dZ.z, 3)-1;
    axes(h(k));
    imagesc(dZ.x, dZ.y,  (dZ.z(:,:,k+1)-dZ.z(:,:,k))/(dZ.t(k+1)-dZ.t(k))); axis xy equal tight; 
    caxis([-4 4]); 
end
end

return

load('/Volumes/insar4/gmap/accum/RACMO2.3/AA//zfirn_20km_10days_all_AA.mat')

K=interp1([-18 0 18], [0 1 0], -18:18);
K=K/sum(K(:));
zfirn.zss=convn(zfirn.zs, reshape(K, [1 1 length(K)]),'same');
zfirn.zss=zfirn.zss./convn(ones(size(zfirn.zs)), reshape(K, [1 1 length(K)]),'same');
[xg1, yg1]=meshgrid(dZ.x, dZ.y);
for k=1:size(dZ.z,3); 
    temp=interpn(zfirn.y, zfirn.x, zfirn.year, zfirn.zss,  yg1, xg1, ones(size(xg1))*dZ.t(k)+2000); 
    dZ.zfirn(:,:,k)=temp; 
end


if 0
    k1=10;
    clf; imagesc(dZ.x, dZ.y, dZ.z(:,:,k1+2)-dZ.z(:,:,k1)),  hold on; for k=1:length(x0); plot(x0(k),'k.','tag', Flist{k}); end;  set(findobj(gca,'type','line'),'buttondownfcn','disp(get(gcbo,''tag''))');
    clf; imagesc(dZ.x, dZ.y, dZ.z(:,:,end-2)-dZ.z(:,:,2)),  hold on; for k=1:length(x0); plot(x0(k),'k.','tag', Flist{k}); end;  set(findobj(gca,'type','line'),'buttondownfcn','disp(get(gcbo,''tag''))');
end

temp=(dZ.season_count>5).*repmat(reshape(dZ.t, [1 1 length(dZ.t)]), [size(dZ.season_count, 1), size(dZ.season_count, 2), 1]);
temp(temp==0)=NaN;
dZ.trange=max(temp, [], 3)-min(temp, [], 3);

dZ.dzdt=(dZ.z(:,:,end-3)-dZ.z(:,:,2))/(dZ.t(end-3)-dZ.t(2));

