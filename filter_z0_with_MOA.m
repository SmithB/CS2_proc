function [MOA_mask, MOA, d_MOA]=filter_z0_with_MOA(D, z0, grids, MOA, d_MOA)

if ~exist('MOA','var') || isempty(MOA)
    MOA=read_geotif_xy('/Volumes/ice1/ben/MOA/2012/moa125_hp1_2013310_2014067_v33.tif', range(grids.z0.ctrs{2})+[-1000 1000], range(grids.z0.ctrs{1})+[-1000 1000]);
    MOA.z=double(MOA.z);
    smooth_scale=diff(grids.z0.ctrs{1}(1:2))/diff(MOA.x(1:2));
    bad_MOA=MOA.z<1200;
    MOA.z(bad_MOA)=NaN;
    MOA.z=conv_corrected(MOA.z, gaussian(-ceil(6*smooth_scale):ceil(6*smooth_scale), 0, smooth_scale), true);
    MOA.z(bad_MOA)=NaN;
    d_MOA=interp2(MOA.x, MOA.y, double(MOA.z)/17000, grids.z0.ctrs{2}(:)', grids.z0.ctrs{1}(:));
end
 good_MOA=isfinite(d_MOA);

[gx, gy]=gradient(z0, grids.z0.ctrs{2}, grids.z0.ctrs{1});
G_slope=double([ones(size(gx(good_MOA))) gx(good_MOA) gy(good_MOA)]);
d_MOA=d_MOA(good_MOA);

tse_good=good_MOA(:);
for kk=1:4
    m_slope=G_slope(tse_good,:)\d_MOA(tse_good);
    r_slope=NaN(size(z0));
    r_slope(good_MOA)=d_MOA-G_slope*m_slope;
    sigma_slope=max(0.01, iqr(r_slope(good_MOA(:) & tse_good(:)))/2);
    tse_good=isfinite(r_slope) & abs(r_slope) < 3*sigma_slope;
end
r_slope(abs(r_slope) > 3*sigma_slope)=NaN;
MOA_mask=interp2(grids.z0.ctrs{2}(:)', grids.z0.ctrs{1}(:), double(isfinite(r_slope)), D.x, D.y)>0.5;

