function D_out=proc_CS2_swath_BD(P, C, Ph, bursts)

dYdPhi=7.5e5*.022084/2/pi/1.17;
phase_gap_size=200/dYdPhi;

Cbad=dilate_mask(C>0.995 | (C>0.900 & P<2e-17), ones(7,1));% eliminate low-signal high-coherence sections (where we suspect that part of the burst has no signal from one RX chain);
C(Cbad)=0;

if 1
    K=gaussian((-16:16)', 0, 5);
    Phc=exp(1i*Ph);
    Cw=conv2(C, K,'same');
    Phs=angle(conv2(Phc.*C, K,'same')./Cw);   
    Cs=conv2(C.*C, K,'same')./Cw;
    %CPh=abs(conv2(Phc.*C, K,'same')).^2./conv2(abs(Phc).^2.*C.^2, K,'same');
    Ph0=Ph;
    Ph=Phs;
    C0=C;
    C=Cs; %.*CPh;
    
end

C_threshold=0.750;

if ~exist('bursts','var') || isempty(bursts)
    bursts=1:size(P,2);
end

if 0
    [good, ind]=ismember(bursts, D_POCA.burst);
    bursts=bursts(good);
    ind=ind(good);
    samp_POCA=D_POCA.samp(ind);
    Ph_POCA=D_POCA.phase(ind);
end

fields={'burst','samp','phase','seg_ind', 'coherence'};
for kF=1:length(fields); 
    seg_dummy_struct.(fields{kF})=[];
end

D_temp=repmat(seg_dummy_struct, 1, length(bursts));
samp=1:size(Ph,1); 
for kB=1:length(bursts);
    k=bursts(kB);
    
    % heal small gaps in the coherence mask
    this_good=conv(double(C(:,k)>C_threshold), ones(5,1),'same')>2;
    
    % unwrap the phase for this trace, just to look for discontinuities 
    % we'll unwrap the phases for the individual segments later on
    temp=unwrap(Ph(:, k));
    temp(~this_good)=NaN;
    
    % find segments with continuous phase
    good_samps=samp(isfinite(diff(temp)) & abs(diff(temp))<phase_gap_size);
    if isempty(good_samps); continue; end
    gap_ind=find_gaps(good_samps, 5);
    clear segs; 
    segs(1,:)=good_samps(gap_ind(1,:)); 
    segs(2,:)=good_samps(gap_ind(2,:));
   
    % get rid of small segments
    segs=segs(:, segs(2,:)-segs(1,:)>20);
   
    if isempty(segs); continue; end
    
    D=repmat(seg_dummy_struct, size(segs,2),1);
    for kSeg=1:size(segs,2);
        % unwrap, starting at the 10th sample (avoids edge problems)
        ind=(segs(1,kSeg):segs(2,kSeg))';
        Ph_seg=NaN(length(ind),1);
        Ph_seg(10:end)=unwrap(Ph(ind(10:end),k));
        Ph_seg(10:-1:1)=unwrap(Ph(ind(10:-1:1),k));
        delta_amb=Ph_seg-Ph(ind,k);
        % check that the phase gradient shows the ground-return point moving away from the POCA
%         if Ph_seg(1) > Ph_POCA(kB) && Ph_seg(end)-Ph_seg(1) < 0
%             Ph_seg=Ph_seg-2*pi;
%         elseif Ph_seg(1) < Ph_POCA(kB) && Ph_seg(end)-Ph_seg(1) > 0
%             Ph_seg=Ph_seg+2*pi;
%         end
        D(kSeg).burst=zeros(size(Ph_seg))+k;
        D(kSeg).samp=ind(:);
        D(kSeg).phase=Ph_seg(:);
        D(kSeg).phase_raw=Ph0(ind,k)+delta_amb;
        D(kSeg).seg_ind=kSeg+zeros(size(Ph_seg));
        D(kSeg).coherence=C(ind, k);   
        D(kSeg).coherence_raw=C0(ind, k);
        D(kSeg).power=P(ind,k);
    end
    
    % stitch together the segments for this burst
    F=fieldnames(D);
    for kF=1:length(F); 
        D_temp(kB).(F{kF})=cat(1, D.(F{kF}));
    end
end

F=fieldnames(D_temp);
for kF=1:length(F);
    D_out.(F{kF})=cat(1, D_temp.(F{kF}));
end

%-------------------------------------------------------------------
function segs=find_gaps(x, tol)

delta=abs(diff(x));
gaps=find(delta>tol);
if isempty(gaps);
    segs=[1; length(x)];
else
    segs=[ 1     gaps(:)'+1 ;
        gaps(:)' length(x)];
    segs=segs(:,diff(segs)>0);
end
