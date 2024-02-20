function D=read_geotif_xy(filename, XR, YR)


% 
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

[x,y, nx, ny, dx, dy]=get_geotif_info(filename);

if length(XR)==3;
    x_skip=max(1, round(XR(2)/dx));
    y_skip=max(1, round(YR(2)/dx));
    XR=XR([1 3]);
    YR=YR([1 3]);
elseif length(XR)==1
   x_skip=max(1, round(XR/dx));
    y_skip=max(1, round(YR/dx));
    XR=x([1 end]);
    YR=y([1 end]);
else
    x_skip=1;
    y_skip=1;
end

cols=find(x>=XR(1) & x<=XR(2)); 
rows=find(y>=min(YR) & y<=max(YR));
if ~any(cols) || ~any(rows); D=struct('x',[],'y',[],'z',[]); return; end
cols=range(cols); cols=[cols(1); x_skip; cols(2)];
rows=range(rows); rows=[rows(1); y_skip; rows(2)];

[D.z, D.x, D.y]=read_geotif(filename, rows, cols);
  

