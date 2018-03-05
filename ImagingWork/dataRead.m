%% Read of Data with adding Missed Data
path = 'G:\MSIBILGISAYAR\UTD\Research\LoFASM\harddiskRecordedData\20170922-20170923_1day_lofasm\csv\AD';
D = dir(path);
addpath(path);
AB = [];
for k = 1:6
    ABtemp = transpose(csvread(D(k+2).name));
    AB = [AB ABtemp];
end

%% IFFT of Each Column
% ABifft = ifft(AB);

%% FFT of Each Column
% AB = fft(ABifft);