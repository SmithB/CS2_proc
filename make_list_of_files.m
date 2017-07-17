 

load /Volumes/insar6/ben/Cryosat/POCA_h5_C/south/AA_helm_v1c_sigma4/bins_for_CS_files.mat
 XR=[-1635000    -1565000]+[-1 1]*1.5e4;
 YR=1.0e+05*[ -2.8500   -2.1500]+[-1 1]*1.5e4;
 x0=mean(XR); dx0=diff(XR);
 y0=mean(YR); dy0=diff(YR);
 
 good=false(length(d),1);
 for k=1:length(X0)
     these=real(X0{k})>x0-dx0/2 & real(X0{k}) <= x0+dx0/2  & imag(X0{k})>y0-dx0/2 & imag(X0{k}) <= y0+dx0/2;
     if any(these)
         good(k)=true;
     end
 end
         
 d=d(good);
 X0=X0(good);
 for k=1:length(d); 
     d(k)=strrep(d(k),'.mat','.DBL');
 end
 
 save list_of_files d
 