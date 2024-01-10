function  varargout=make_AA_h5(varargin)

if nargout>0
    [varargout{1:nargout}]=feval(varargin{:});
else
    feval(varargin{:});
end


function queue_list=make_queue(in_dir, out_dir, dx0, mode);

%in_dir='slope_picks/south/AA_v2s/'; out_dir='/Volumes/ice1/ben/Cryosat/SW_h5/south/AA_v2s/';make_AA_sw_h5(in_dir, out_dir);
%in_dir='slope_picks_C/south/AA_v1b_sigma4/'; out_dir='/Volumes/insar6/ben/Cryosat/SW_h5_C/south/AA_v1b_sigma4/' ;make_AA_sw_h5(in_dir, out_dir);


% first: collect the x0 variables from the input files, make the queue for
% the indexing
%xr=[   -700000   900000];
%yr=[   -3300000  -700000 ];

 
switch mode
    case 'sw'
        [~, d]=unix(['ls ', in_dir,' | grep _sw.mat']);
        prefix='sw_';
    case 'POCA'
        [~, d]=unix(['ls ', in_dir,' | grep -P "CS.*00..mat"']);
        prefix='POCA_';
end

d=strsplit(deblank(d));
for k=1:length(d);
    d{k}=strrep(d{k}, in_dir,'');
    temp=load([in_dir, '/', d{k}],'X0');
    if ~isempty(temp);
        X0{k}=temp.X0;
    end
end

if ~exist(out_dir,'dir');
    mkdir(out_dir);
end
uBins=unique(round_to(cat(1, X0{:}), dx0));
fid=fopen([out_dir,'/bin_list.txt'],'w');
for k=1:length(uBins);
    fprintf(fid,'%s%dE_%dN.h5\n', prefix, real(uBins(k))/1000, imag(uBins(k))/1000);
end
fclose(fid);


save([out_dir,'/bins_for_CS_files.mat'], 'X0', 'd', 'dx0', 'in_dir');

fid_list=fopen([out_dir,'/bin_list.txt'],'r');

while ~feof(fid_list)
    fname=deblank(fgetl(fid_list));
    queue_list{k}=sprintf( 'make_AA_h5(''%s'',''%s'', ''%s'');', in_dir, out_dir, fname);
end
fclose(fid_list);

if nargout==0
    fid_queue=fopen(sprintf('queue_for_AA%s.txt',mode),'w');
    for k=1:length(queue_list)
        fprintf(fid_queue, '%s\n', queue_list{k});
    end
    fclose(fid_queue);
end


function run_bin(in_dir, out_file, mode)

% I am a worker. 
if ~isempty(in_dir) && exist('out_file','var'); 
    fprintf(1,'make_AA_h5: working on %s\n', out_file);
    load([out_dir,'/bins_for_CS_files.mat'])
    temp=regexp(out_file,'([-+]*\d+)E_([-+]*\d+)N','tokens');
    x0=str2double(temp{1}{1})*1000;
    y0=str2double(temp{1}{2})*1000;
    
    clear D0
    for k=1:length(X0);
        these=real(X0{k})>x0-dx0/2 & real(X0{k}) <= x0+dx0/2  & imag(X0{k})>y0-dx0/2 & imag(X0{k}) <= y0+dx0/2;
        if ~any(these);
            continue
        end
        
        try
            L=load([in_dir,'/',d{k}]);
            switch mode
                case 'sw'
                    if isfield(L,'Dsw') && ~isempty(L.Dsw);
                        D0(k)=L.Dsw;
                        temp_flag=D0(k).POCA_flag;
                        D0(k).POCA_flag=false(size(D0(k).burst));
                        D0(k).POCA_flag(1:length(temp_flag))=temp_flag;
                    end
                case 'POCA'
                    
            end
        catch
            fprintf(1,'Failure for %s\n', [in_dir,'/',d{k}]);
        end
    end
    if ~exist('D0','var'); return; end
    f=fieldnames(D0);
    for kf=1:length(f)
        D.(f{kf})=cat(1, D0.(f{kf}));
        [D0.(f{kf})]=deal([]);  % try to control memory use
    end
    these=real(D.x)>x0-dx0/2 & real(D.x) <= x0+dx0/2  & imag(D.x)>y0-dx0/2 & imag(D.x) <= y0+dx0/2 & isfinite(D.x);
    if ~any(these); return; end
    D=index_struct(D, these);
    
    % check this
    [M, INDEX]=index_point_data_h5('build_index_and_sort', D, [], out_dir, 1e4, 10.2);
    
    index_point_data_h5('write_h5',[out_dir,'/',out_file], M, INDEX, true);
    save([out_dir,'/', strrep(out_file,'.h5','_index.mat')], 'INDEX');
end


function make_index(in_dir, out_dir);
% If the input directory argument is empty, we've made the subset h5 files, now 
% make the master indeaddx
if isempty(in_dir) && exist(out_dir,'dir');
    [~, h5_files]=unix(['ls ', out_dir,'/*.h5']);
    h5_files=strsplit(deblank(h5_files));
    for k=1:length(h5_files);
        try
            II=h5info(h5_files{k},'/INDEX');
            disp(II.Filename);
        catch
            fprintf(1,'failure for %s\n', h5_files{k});
            delete(h5_files{k})
        end
    end
    
    INDEX_list=index_point_data_h5('collect_index', {out_dir}, '*.h5',  {'./'});
    
    master_index=index_point_data_h5('make_master_index', INDEX_list, [out_dir,'/master_index_h5.mat'],  out_dir);
end