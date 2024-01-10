function D=select_swath_ambiguity(D_in, DEM, out_fieldnames)

DOPLOT=false;
if ~exist('out_fieldnames','var')
    out_fieldnames={'xPS','h','ambiguity', 'burst','seg_ind', 'phase','power','coherence','time','samp'};
end
if ~isfield(D_in,'ambiguity')
    Dsw.ambiguity=repmat([-1 0 1], [size(D_in.h,1),1]); Dsw.ambiguity=Dsw.ambiguity(:);
end
dYdPhi=7.5e5*.022084/2/pi/1.17;

if ~isfield(D_in,'seg_ind')
    D_in.seg_ind=D_in.ret_count;
end

% replicate fields that need it (to cover ambiguities) and flatten
for this=[out_fieldnames,{'seg_ind'}]   
    if strcmp(this{1},'ambiguity'); continue; end
    temp=D_in.(this{1});
    if size(temp,2)==1
        temp=repmat(temp, [1,3]);
    end
    Dsw.(this{1})=temp(:);
end

Dsw.DEM=NaN(size(Dsw.h));
if isfield(DEM,'t0')
    Dsw.DEM=interpn(DEM.y, DEM.x, [DEM.t0, DEM.t0(end)+median(diff(DEM.t0))], cat(3, DEM.z, DEM.z(:,:,end)), imag(Dsw.xPS), real(Dsw.xPS), Dsw.time, '*linear')+ ...
        interpn(DEM.y, DEM.x, DEM.z0,imag(Dsw.xPS), real(Dsw.xPS), '*linear') ;
else
    if min(size(DEM.z)) > 2
        Dsw.DEM=interp2(DEM.x, DEM.y, DEM.z, real(Dsw.xPS), imag(Dsw.xPS),'*linear');
    end
end
if isfield(Dsw,'geoid') & any ~isfinite(Dsw.DEM);
    Dsw.DEM(~isfinite(Dsw.DEM))=Dsw.geoid(~isfinite(Dsw.DEM));
end
Dsw.r=Dsw.h-Dsw.DEM;

amb_list=-1:1;
 
if 0   
    XR=[-1.635 -1.565]*1e6;
    YR=[-3.5 -2.35]*1e5;
    % to select a file:
    if 0
        mat_wc='/Volumes/ice1/ben/Cryosat/slope_picks/south/Amundsen_v3s/*_sw.mat';
        [~, out]=unix(['ls ', mat_wc]);
        out=strsplit(deblank(out));
        InPolyCount=zeros(size(out));
        XR1=XR+[-1e4 1e4]; YR1=YR+[-1e4 1e4];
        for k=1:length(out)
            load(out{k},'X0');
            InPolyCount(k)=sum(real(X0) >= XR1(1) & real(X0) <= XR1(2) & imag(X0) >= YR1(1) & imag(X0) <= YR1(2));
        end      
    end
    these=real(Dsw.xPS)> XR(1) & real(Dsw.xPS) < XR(2) &imag(Dsw.xPS)> YR(1) & imag(Dsw.xPS) < YR(2) ;
    uB=unique(Dsw.burst(these));
    ff={'xv','zv','h_sc','xps_sc','yps_sc'};
    for k=1:length(ff); Dsw=rmfield(Dsw, ff{k}); end 
    Dsw=index_struct(Dsw, ismember(Dsw.burst, uB));
end


BS=Dsw.burst(:)+0.0001*Dsw.seg_ind(:);
% loop over distinct bursts and segment indices
els=bin_by( BS, unique(BS));
D_out=repmat(struct('count', [], 'dh_spline_dx', [], 'block_h_spread', [],'index', []), [1,length(els)] ); 

for kB=1:length(els)
    D1=index_struct(Dsw, els{kB},{'phase', 'h','r','ambiguity','power','coherence'});
    D1.index=els{kB};
    if all(D1.phase==0); continue; end
    R=NaN(3,1);
    spread=NaN(3,1);
    N=NaN(3,1);
    for kA=1:3
        these=D1.ambiguity==amb_list(kA) & isfinite(D1.r);
        if any(these)
            R(kA)=mean(abs(D1.r(these)));
            spread(kA)=mean(abs(D1.r(these)-median(D1.r(these))));
            N(kA)=sum(these);
        end
    end
    if all(isnan(R)); continue; end
    best_R=find(isfinite(R) & R==min(R(isfinite(R))),1,'first');
    best_spread=find(isfinite(spread) & spread==min(spread(isfinite(spread))),1,'first');
    if all(N(best_R)>3 & N(best_spread)>3)
        % swath case
        if best_R ~=best_spread  
            % calculate the fractional improvement in spread and range from
            % selecting the smallest, rather than the next-smallest
            % ambiguity
            frac_spread=(spread-min(spread))./max(1,min(spread));
            frac_R=(R-min(R))./max(2,min(R));
             
            if frac_spread(best_R) < 0.5           
                best=best_R;
            elseif frac_R(best_spread) <0.1
                % there's only a 10% range improvement from selecting this
                % ambiguity
                best=best_spread;
            else
                % we can't tell which ambiguity is best, punt
                continue
            end
        else
            % spread and R tell the same story
            best=best_R;
        end
    else
        % single-point case
        best=best_R;
    end
    
    amb_best=amb_list(best);
    R_sorted=sort(R);
    D2=index_struct(D1,D1.ambiguity==amb_best & D1.coherence > 775 );
    D2.delta_ambiguity=zeros(size(D2.h))+(R_sorted(2)-R_sorted(1));
    if isempty(D2.h); continue; end
    x_phase=(D2.phase+2*pi*D2.ambiguity)*dYdPhi;
    if length(D2.h) > 1        
        if length(D2.h) > 3
            x_node=linspace(min(x_phase)-50, max(x_phase)+50, max(2,ceil(diff(range(x_phase))/400)));
            if length(x_node) == 2
                G_spline=[x_node(2)-x_phase(:) x_phase(:)-x_node(1)]/diff(x_node);
            else
                G_spline=fd_spline_fit('build_G', x_phase, x_node, 0);
            end
            h_node=[G_spline; 0.0001*(eye(size(G_spline,2))-(1/size(G_spline,2)))]\[D2.h; zeros(size(G_spline,2),1)];
            block_r=D2.h-G_spline*h_node;
            
            dh_spline_dx=diff(h_node)./diff(x_node(:));
            
        else
            block_r=D2.h;
            dh_spline_dx=diff(D2.h([1 end]))./diff(x_phase([1 end]));
            x_node=x_phase([1 end]);
        end
        [~, ~, bin_percentiles, count]=block_percentile_filter(double(x_phase), double(block_r), [.16 .84], 400, [0 200]);
        ind=blockmedian(double(x_phase), double(block_r), 400, [0 200]);
        [~, temp]=sort(x_phase(ind)); ind=ind(temp); bin_percentiles=bin_percentiles(temp,:); count=count(temp);
        D_out(kB).count=count;
        D_out(kB).block_h_spread=bin_percentiles(:,2)-bin_percentiles(:,1);
        if length(x_node)>2
            D_out(kB).dh_spline_dx=interp1((x_node(1:end-1)+x_node(2:end))/2, dh_spline_dx, x_phase(ind),'linear','extrap');
        else
            D_out(kB).dh_spline_dx=dh_spline_dx(1)*ones(size(ind)); 
        end
        
        
        if DOPLOT
            figure(5); plot(x_phase, D2.h,'--', x_phase, G_spline*h_node,'r.',...
                x_node, h_node,'g*',...
                x_phase(ind), D2.h(ind), ...
                x_phase(ind), D2.h(ind)-D_out(kB).block_h_spread/2, ...
                x_phase(ind), D2.h(ind)+D_out(kB).block_h_spread/2); 
            set(findobj(gca,'type','line','linestyle','-'),'marker','.')
            disp('next...');
        end
    else
        ind=1;
        D_out(kB).count=1;
        D_out(kB).block_h_spread=0;
        D_out(kB).dh_spline_dx=NaN;
    end
    D_out(kB).index=D2.index(ind);
end

if ~exist('D_out','var')
    D=[];
    return; 
end
% cat together D_out;
D_out_fieldnames=fieldnames(D_out);
for kf=1:length(D_out_fieldnames)
    D.(D_out_fieldnames{kf})=cat(1, D_out.(D_out_fieldnames{kf}));
end

% copy data from Dsw to D_out
for kf=1:length(out_fieldnames)
    D.(out_fieldnames{kf})=Dsw.(out_fieldnames{kf})(D.index);
end

D=rmfield(D, 'index');

