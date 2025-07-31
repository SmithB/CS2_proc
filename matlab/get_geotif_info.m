function [x,y, nx, ny, dx, dy]=get_geotif_info(filename);


% unix_cmd=['gdalinfo ', filename,'| perl -e ''', ...
%     'while($line =<>) {if ($line =~ /Origin\s*=\s*\(\s*(\S*\d+.\d+),(\S*\d+.\d+)/){print "$1 $2 "}', ...
%     'elsif (($line =~ /Pixel Size/) & ($line=~ /\((\S*\d+.\d+),(\S*\d+.\d+)/)) { print "$1 $2\n"}', ...
%     'elsif ($line =~ /ize\s*is\s*(\d+),\s*(\d+)/) {print "$1 $2 "}}'''];
% 
% [s,info]=unix(unix_cmd);
% 
% info=str2num(info);
% nx=info(1); ny=info(2);
% ulx=info(3); uly=info(4);
% dx=info(5); dy=info(6);
% x=ulx+(0:nx-1)*dx;
% y=uly+(0:ny-1)*dy;


% Yushin Ahn's code to read tags:
Tinfo        = imfinfo(filename);
Tinfo=Tinfo(1);
info.samples = Tinfo.Width;
info.lines   = Tinfo.Height;
info.imsize  = Tinfo.Offset;
info.bands   = Tinfo.SamplesPerPixel;


sub = [1, info.samples, 1, info.lines];
%data_type = Tinfo.BitDepth/8;
data_type = Tinfo.BitsPerSample(1)/8;
switch data_type
    case {1}
        format = 'uint8';
    case {2}
        format = 'int16';
    case{3}
        format = 'int32';
    case {4}
        format = 'single';
end

try
    info.map_info.dx = Tinfo.ModelPixelScaleTag(1);
    info.map_info.dy = Tinfo.ModelPixelScaleTag(2);
    info.map_info.mapx = Tinfo.ModelTiepointTag(4);
    info.map_info.mapy = Tinfo.ModelTiepointTag(5);
catch
    info.map_info.dx=Tinfo.ModelTransformationTag(1);
    info.map_info.dy=Tinfo.ModelTransformationTag(6);
    info.map_info.mapx = Tinfo.ModelTransformationTag(4);
    miny=Tinfo.ModelTransformationTag(8);
    info.map_info.mapy = miny+(info.samples-1)*info.map_info.dy;
end
%info.map_info.projection_name = Tinfo.GeoAsciiParamsTag;
%info.map_info.projection_info = Tinfo.GeoDoubleParamsTag;

minx = info.map_info.mapx;
maxy = info.map_info.mapy;
maxx = minx + (info.samples-1).*info.map_info.dx;
miny = maxy - (info.lines-1  ).*info.map_info.dy;


xm = info.map_info.mapx;
ym = info.map_info.mapy;
% N.B.  -- modified 12/26/12 to report values at pixel centers
x = xm + ((0.5:info.samples-0.5).*info.map_info.dx);
y = ym - ((0.5:info.lines  -0.5).*info.map_info.dy);

nx=info.samples; 
ny=info.lines; 
dx=info.map_info.dx;
dy=info.map_info.dy;