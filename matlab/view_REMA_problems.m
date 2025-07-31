[~, files]=unix('ls Problems/1/DEMs');
files=strsplit(deblank(files));
clear I I1 temp;
for k=1:length(files); I(k)=read_geotif_xy(['Problems/1/DEMs/', files{k}], range(F.ll(:, 1))+[-2e4 2e4], range(F.ll(:,2))+[-2e4 2e4]); end

for k=1:length(I); 
    [xg, yg]=meshgrid(I(k).x(1:20:end), I(k).y(1:20:end)); 
    zg=I(k).z(1:20:end, 1:20:end); 
    temp(k)=struct('x', xg(isfinite(zg)),'y', yg(isfinite(zg)),'z', zg(isfinite(zg))); 
end
temp=flatten_struct(temp);
x0=mean(temp.x); y0=mean(temp.y);sx=std(temp.x); sy=std(temp.y);

temp.x=(temp.x-x0)/sx; temp.y=(temp.y-y0)/sy;
G=[ones(size(temp.x)), temp.x temp.y temp.x.*temp.y temp.x.^2 temp.y.^2];
m=G\temp.z;

figure(1); clf; hold on;
I1=I;
for k=1:length(I)
    [xg,yg]=meshgrid((I(k).x-x0)/sx, (I(k).y-y0)/sy);
    xg=xg(:); yg=yg(:);
    I1(k).z=I1(k).z-reshape([ones(size(xg)), xg yg xg.*yg xg.^2 yg.^2]*m, size(I1(k).z));
    hs(k)=surf(I1(k).x(1:5:end), I1(k).y(1:5:end), I1(k).z(1:5:end, 1:5:end),'tag', files{k}); shading interp;
    
end

set(findobj(gca,'type','surface'),'buttondownfcn','disp(get(gcbo,''tag''))')
