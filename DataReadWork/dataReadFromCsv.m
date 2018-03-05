addpath(genpath('C:\D\UTD\Research\LoFASM\harddiskRecordedData'));

%% Read of Data with adding Missed Data
path = 'C:\D\UTD\Research\LoFASM\harddiskRecordedData\LoFASM1\csv\DD';
D = dir(path);
addpath(path);
R = [];
for k = 1:105
    Rtemp = transpose(csvread(D(k+2).name));
    R = [R Rtemp];
end

%% Data Initial Analysis
% addpath(genpath('C:\D\UTD\Research\LoFASM\ImagingWork'));
% RwF = whiteningApproach1(R);
% 
% RlTerm = longTermAverage(RwF,100);
% mesh(mag2db(abs(RlTerm)))
% 
% RwFifft = ifft(RlTerm);
% figure;mesh(abs(RwFifft))