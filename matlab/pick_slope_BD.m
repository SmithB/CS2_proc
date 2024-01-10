function D_POCA=pick_slope(P, C, Ph, bursts, noise_dB)

bad=(C>2); C(bad)=NaN;

% P is the L1B power matrix, reshaped to N_samples x N_traces
% C is the L1B coherence matrix, same shape as P
% Ph is the phase matrix, same shape as P, in units of radians (scale raw
% values by 1/1e6);
 
%S0=max( max(P)/1e3, std(P(1:20,:)));
if exist('noise_dB','var')
    S0=10.^(noise_dB(:)/10)';
else
    S0=5e-17*ones(size(P(1,:)));  % there is no evidence that this should vary.  It's probably a good idea to look for a better representative value...
end

% NOTE FROM BEN: Some traces have (effectively) a zero noise level.
% for these, replace the noise level with the median noise level from all
% traces.
S0(S0<1e-20)=median(S0(S0>1e-20)); %%%

% width of the smoothing kernel
% This was set to 4 for the smooth pick set.  
% a value of 2 can find some sharper returns, may be less biased!
sigma_k=4;


% approximate scaling from phase difference to cross-track distance : H
% lambda / (2 pi B), ~2200m/rad
dYdPhi=7.5e5*.022084/2/pi;

% build a smoothing kernel
K_smooth=gaussian(-16:16, 0, sigma_k);
% build a slope-smoothing kernel
K_smooth_slope=conv(K_smooth, [1 0 -1]/2,'same');
% calculate the smooth slope
dPdt=conv2(P, K_smooth_slope(:), 'same');

% make a mask for the data in which slopes are significantly different from
% zero and where the coherence is high enough (> 750)
mask=dPdt > repmat(6*S0/sqrt(sigma_k), [size(P,1),1]) & C>0.750;

% kill off small islands in the mask
mask=conv2(double(conv2(double(mask), ones(sigma_k+1,1),'same')>sigma_k/2), ones(sigma_k,1), 'same')>0;

Smax_ind=deal(NaN(size(P,2),1));

delta_ind=(-1:1)';
G=[ones(size(delta_ind)), delta_ind, delta_ind.^2];
Ginv=(G'*G)\G';

if ~exist('bursts','var') || isempty(bursts)
    bursts=1:size(P,2);
end

[P_est, phase_est, C_est, P_int, N_FP_bins, dPdt_est, P_pre, burst_num, ret_count]=deal(NaN(size(P,2),1));

fields={'power','coherence','burst','samp','phase','P_int','N_fp_bins','dPdt_est','power_before_pick', 'ret_count'};
for kF=1:length(fields)
    D_POCA.(fields{kF})=[];
end

count=0;
for ind=1:length(bursts)
    k=bursts(ind);
    
    % find the maxima in dPdt
    maxima=1+find(dPdt(2:end-1,k)>dPdt( 1:end-2,k) & dPdt(2:end-1,k)> dPdt(3:end,k));
    
    if isempty(maxima) ||~any(mask(maxima,k))
        continue
    end
    
    % select the maxima that are in the mask
    maxima=maxima(mask(maxima,k));
    % smooth the power so that we can interpolate a robust value from it
    P_smooth=conv(P(:,k), K_smooth(:),'same');
    C_smooth=conv(C(:,k), K_smooth(:),'same');
    
    P_maxima=1+find(P_smooth(2:end-1)>P_smooth( 1:end-2) & P_smooth(2:end-1)> P_smooth(3:end));
    P_maxima=P_maxima(P_smooth(P_maxima) > 0.5*sqrt(mean(P_smooth).^2));
    
    maxima=maxima(dPdt(maxima,k) > 0.5*(P_smooth(maxima)/2/sigma_k));
    P_envelope=interp1([0; P_maxima(:)], P_smooth([P_maxima(1); P_maxima(:)]), maxima);
    maxima=maxima(P_smooth(maxima)>0.2*P_envelope);
    
    if isempty(maxima); continue; end
    clear reg1
    for k_max=1:length(maxima)
        this_max_ind=maxima(k_max);
        % select a range of samples around this maximum
        i0=find(mask(1:this_max_ind,k)==0, 1, 'last');
        if isempty(i0)
            i0=1;
        else
            i0=i0+1;
        end
        
        
        i1=find(~mask(i0:end,k), 1, 'first');
        if isempty(i1) || i1 < 5
            continue
        end
        
        % we've found something...
        count=count+1;
        
        %report this maximum location
        Smax_ind(count)=this_max_ind;
        burst_num(count)=k;
        
        
        % do a parabolic refinement to get a refined estimate of the
        % maximum-slope sample
        if Smax_ind(count) < i1
            m=Ginv*dPdt(Smax_ind(count)+delta_ind,k);
            Smax_ind(k)=Smax_ind(count)-m(2)/2/m(3);
        end
        % pull out the phase for the region around Smax_ind
        reg1{k_max} = (floor(Smax_ind(count))+(-6*sigma_k:6*sigma_k));
        reg1{k_max} = reg1{k_max}(reg1{k_max}>0 & reg1{k_max} < size(P,1));
        Ph1=unwrap(Ph(reg1{k_max},k));
        
        % since ambiguities can pile up in the unwrapping process, want to
        % anchor the unwrapping to the phase of the first good bin (i.e. the
        % first bin in the mask for this trace)
        delta_phi_unwrap=Ph1-Ph(reg1{k_max}, k);
        first_good_bin=find(mask(reg1{k_max},k), 1, 'first');
        if ~isempty(first_good_bin)
            Ph1=Ph1-delta_phi_unwrap(first_good_bin);
        end
        
        % smooth the phase
        Ph1=conv(Ph1, K_smooth,'same')./conv(ones(size(Ph1)), K_smooth,'same');
        phase_est(count)=interp1(reg1{k_max}, Ph1, Smax_ind(count));
        
        P_est(count)=interp1(1:length(P_smooth), P_smooth, Smax_ind(count));
        
        dPdt_est(count)=interp1(reg1{k_max}, dPdt(reg1{k_max},k), Smax_ind(count));
        
        first_samp=Smax_ind(count)-2*P_est(count)/dPdt_est(count);
        if first_samp > 1
            P_pre(count)=mean(P(1:floor(first_samp),k));
        else
            P_pre(count)=NaN;
        end
               
        % find elements inside a 300 m-radius footprint around the return location and
        % add up their power.  The maximum range window for the FP is
        % equivalent to 4x the inverse bandwidth (8 bins on either side)
        
        els_in_fp=dYdPhi*abs(Ph1-phase_est(count)) < 300 & abs(reg1{k_max}(:) - Smax_ind(count)) < 8;
        P_int(count)=sum(P(els_in_fp,k));
        N_FP_bins(count)=sum(els_in_fp);
        
        C_est(count)=interp1(1:length(C_smooth), C_smooth, Smax_ind(count));    
        ret_count(count)=k_max;
    end
end    

D_POCA.power=P_est(:);
D_POCA.coherence=C_est(:);
D_POCA.burst=burst_num(:);
D_POCA.samp=Smax_ind(:);
D_POCA.phase=phase_est(:);
D_POCA.P_int=P_int(:);
D_POCA.N_fp_bins=N_FP_bins(:);
D_POCA.dPdt_est=dPdt_est(:);
D_POCA.power_before_pick=P_pre(:);
D_POCA.ret_count=ret_count(:);

D_POCA=index_struct(D_POCA, isfinite(D_POCA.burst+D_POCA.samp));


 