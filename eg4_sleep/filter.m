cd d:\sleep\capslpdb
db_list=physionetdb('capslpdb');

pkg load signal
file2= regexprep(db_list{1,2},'.edf','.mat');
load(file2);

clf
ylimrg= [-20 20];
t1_all= [[0;30;60],[1860;1890;1920], [2310;2340;2370], [4650;4680;4710], [5040;5070;5100], [7590;7620;7650]];
n_type= size(t1_all)(2);
n_sample= size(t1_all)(1);
len_sample= 30; # seconds window
rate_sample= 256; # Herz
fc = [0 20]; # cutoff frequency
fs = 256; # sampling frequency
[b,a] = butter(1,fc/(fs/2)); # https://www.mathworks.com/help/signal/ref/butter.html
n = len_sample*rate_sample;        # frequency spectrum  
t_power = (0:n-1)*(fs/n);     # frequency range
ncol= 3;

for f= 1:n_type
  figure(f)
  data= zeros(len_sample*rate_sample,n_sample);
  for s= 1:n_sample
##  raw EEG
    subplot(n_sample,ncol,s*ncol-2)
    low= t1_all(s,f);
    high= low+len_sample;
    x1= temp((low*rate_sample+1):(high*rate_sample),2);
    t1= temp((low*rate_sample+1):(high*rate_sample),1);
    plot(t1,x1)
    xlabel('Time')
    ylabel('EEG voltage')
    ylim(ylimrg)

    data(:,s)= x1;
   end
   filtered = filter(b,a,data);

##   plotting filtered
   for s= 1:n_sample
    subplot(n_sample,ncol,s*ncol-1)
    plot(filtered(:,s))
    xlabel('Time')
    ylabel('Filtered')
    ylim(ylimrg)
    
    ##  frequency spectrum
    subplot(n_sample,ncol,s*ncol)
    y = fft(filtered(:,s));
    power = abs(y).^2/n;    % power of the DFT
    plot(t_power,power)
    xlabel('Frequency')
    ylabel('Power')
    xlim([0 20])

   end
end



####https://www.mathworks.com/help/matlab/math/basic-spectral-analysis.html
##fs = 256;                                % sample frequency (Hz)
##t = 0:1/fs:30-1/fs;                      % 10 second span time vector
##x = x1;         
##y = fft(x);
##n = length(x);          % number of samples
##f = (0:n-1)*(fs/n);     % frequency range
##power = abs(y).^2/n;    % power of the DFT
##figure;
##plot(f,power)
##xlabel('Frequency')
##ylabel('Power')
##xlim([0 20])


##subplot(n_plots,2,1)
##low= 0;
##high= low+30;
##x1= temp.temp((low*256+1):(high*256),2);
##t1= temp.temp((low*256+1):(high*256),1);
##plot(t1,x1)
##ylim(ylimrg)
##
##subplot(2,2,3)
##low= 1860;
##high= low+30;
##x2= temp.temp((low*256+1):(high*256),2);
##t2= temp.temp((low*256+1):(high*256),1);
##plot(t2,x2)
##ylim(ylimrg)
##
##data= [x1, x2];
##
##fc = [2 7]; # cutoff frequency
##fs = 256; # sampling frequency
##[b,a] = butter(1,fc/(fs/2)); # https://www.mathworks.com/help/signal/ref/butter.html
##filtered = filter(b,a,data);
##
##subplot(2,2,2)
##plot(t1,filtered(:,1))
##ylim(ylimrg)
##subplot(2,2,4)
##plot(t2,filtered(:,2))
##ylim(ylimrg)
##
##
##
##
##
## sf = 800; sf2 = sf/2;
## data=[[1;zeros(sf-1,1)],sinetone(25,sf,1,1),sinetone(50,sf,1,1),sinetone(100,sf,1,1)];
## [b,a]=butter ( 1, 5 / sf2 );
## filtered = filter(b,a,data);
##
## clf
## subplot ( columns ( filtered ), 1, 1)
## plot(filtered(:,1),";Impulse response;")
## subplot ( columns ( filtered ), 1, 2 )
## plot(filtered(:,2),";25Hz response;")
## subplot ( columns ( filtered ), 1, 3 )
## plot(filtered(:,3),";50Hz response;")
## subplot ( columns ( filtered ), 1, 4 )
## plot(filtered(:,4),";100Hz response;")
## figure(2)
## plot(data(:,3))
## plot(data(:,4))