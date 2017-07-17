function [count_swath, count_POCA]=count_swath_poca(D, grids)    

D0=D;
count_swath=zeros(grids.z0.dims(1:2));

D=index_struct(D0,D0.swath==1);
bin=round((D.x-grids.z0.ctrs{2}(1))/diff(grids.z0.ctrs{2}(1:2)))+...
    1i*round((D.y-grids.z0.ctrs{1}(1))/diff(grids.z0.ctrs{1}(1:2)));

if ~isempty(D.x);
    
    s=sort(bin);
    [ub, ia, ib]=unique(s,'first');
    [ub, ia1, ib1]=unique(s,'last');
    count_swath(sub2ind(size(count_swath), imag(ub)+1, real(ub)+1))=ia1-ia+1;
end

count_POCA=zeros(grids.z0.dims(1:2));

D=index_struct(D0,D0.swath==0);
bin=round((D.x-grids.z0.ctrs{2}(1))/diff(grids.z0.ctrs{2}(1:2)))+...
    1i*round((D.y-grids.z0.ctrs{1}(1))/diff(grids.z0.ctrs{1}(1:2)));

if ~isempty(D.x);
    
    s=sort(bin);
    [ub, ia, ib]=unique(s,'first');
    [ub, ia1, ib1]=unique(s,'last');
    count_POCA(sub2ind(size(count_POCA), imag(ub)+1, real(ub)+1))=ia1-ia+1;
end