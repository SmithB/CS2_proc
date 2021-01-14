function F=read_nc_field(file, field, start, count, stride)

if nargin > 2
    F=ncread(file, field, start, count, stride);
else
    F=ncread(file, field);
end
try
    bad=F==ncreadatt(file, field, '_FillValue');
catch
    bad=false(size(F));
end
F=double(F);
F(bad)=NaN;
if ~(isa(F,'double') || isa(F, 'single'))
    try
        F=F*ncreadatt(file, field, 'scale_factor');
    catch
    end
end