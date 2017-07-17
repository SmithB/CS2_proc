function roll_corrected=read_roll_correction(burst_time)

roll_dir='/Volumes/ice1/ben/Cryosat/RollBiasCorrection/';
roll_index_file=[roll_dir,'/Roll_index.h5'];
TR=range(burst_time(isfinite(burst_time(:)) & burst_time(:)~=datenum('jan 1 2000')));

IND=h5read(roll_index_file,'/Index');
IND1=find(IND(1,:) < TR(1), 1, 'last');
IND2=find(IND(2,:) > TR(2), 1, 'first');

if isempty(IND1) || isempty(IND2)
    roll_corrected=NaN(size(burst_time)); 
    return
end

file_numbers=min([IND1, IND2]):max([IND1, IND2]);

for k=1:length(file_numbers)
    temp=h5readatt(roll_index_file,'/Index', sprintf('file_%d', file_numbers(k)));
    thefile=[roll_dir,'/h5/',temp{1}];
    S(k)=struct('time', h5read(thefile,'/time'),'roll', h5read(thefile,'/roll'));
end

S=struct('time', cat(1, S.time),'roll', cat(1, S.roll));
[~, ind]=unique(S.time);
S.time=S.time(ind); S.roll=S.roll(ind);
roll_corrected=interp1(S.time, S.roll, burst_time);
