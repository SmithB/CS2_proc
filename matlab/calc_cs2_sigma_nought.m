function Pcorr=calc_cs2_sigma_nought(D)

lambda = 0.022084; %(m)
baseline = 1.1676; %(m)
alpha = .973*lambda*(D.phase+2*pi*D.ambiguity)/(2*pi*baseline); %alpha = inferred angle (rad)  

alpha_3db=1.1992*pi/180;
sigma_beam=alpha_3db*sqrt(-1/(2*log(10^-.3)));

G=exp(-(alpha.^2/2/sigma_beam.^2));

PR=D.range_surf.^2;

Pcorr=D.power./G.^4.*PR;