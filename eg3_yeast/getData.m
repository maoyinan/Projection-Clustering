cd d:\sleep\capslpdb


##file= 'sleep-edfx/sleep-cassette/SC4001E0-PSG.edf'
##x=wfdbdesc(file);
##x.Description
##[sig1, Fs1, tm1] = rdsamp('sleep-edfx/sleep-cassette/SC4001E0-PSG.edf', 1);
##[sig2, Fs2, tm2] = rdsamp('sleep-edfx/sleep-cassette/SC4001EC-Hypnogram.edf');
##figure(1)
##plot(tm1,sig1);
##figure(2)
##plot(tm2,sig2);
##hypnograph doesn't make sense!!!

##try another database CAP 
db_list=physionetdb('capslpdb');
save('db_list.mat', '-mat', 'db_list');

% run ScoringReader.m
% hyp___________containing the hypnogram evaluated for each 30 s epochs
% time_tot______containing the starting time in seconds of CAP phases A
% duration______containing the duration of each phase A in seconds
% type_ar_______containing the type of phase A (A1, A2 or A3)


for iii = setdiff(1:size(db_list)(:,2), [51, 108])
  file= regexprep(db_list{1,iii},'.edf','.txt');
    urlFile= strcat('http://physionet.org/physiobank/database/capslpdb/', file);
    urlwrite(urlFile,file)
    ScoringReader

##    figure(1)
##    plot(hyp(:,2),hyp(:,1));
##    figure(2)
##    plot(time_tot,type_ar);

    file2= strcat('hyp_',regexprep(file,'.txt','.mat'));
    save(file2,'hyp','-v7');
##   load(file2);
end



##read all subjects Fp2-F4 signals
for i = 1:1#size(db_list)(:,2)
    file= strcat('capslpdb/', db_list{1,i});
    temp=wfdbdesc(file);
    ch=find(strcmp({temp.Description},'Fp2-F4')==1);
    fprintf('Reading %s, Fp2-F4 in channel %.0f\n', file,ch);
    [sig, Fs, tm] = rdsamp(file, ch);
##    plot(tm,sig);
    
    temp= [tm(:), sig(:)];
    file2= regexprep(db_list{1,i},'.edf','.mat');
    save(file2,'temp','-v7');
##    temp= load(file2);
end

