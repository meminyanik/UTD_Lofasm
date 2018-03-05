tStart = 1;
fStart = 550; % 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz
fStop = 550; % 1024 = 100 MHz, 820 = 80 MHz, 417 = 40 MHz

SimulationMode = false;

% load('tgExact_LST_50Sample.mat')
load('tgExact_LST.mat')

% load('peakPhaseInFrequencyCygAnoWhitening');
% peakPhaseInFrequency = unwrap(peakPhaseInFrequency);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Filtered Data
%-------------------------------------------------------------------------%
% AAlT = reshape(RFiltered(1,1,:,:),1024,504);
% ABlT = reshape(RFiltered(2,1,:,:),1024,504);
% BBlT = reshape(RFiltered(2,2,:,:),1024,504);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%-------------------------------------------------------------------------%
if SimulationMode
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Simulated Data
    %-------------------------------------------------------------------------%
%     load('AASim');
%     load('ABSim');
%     load('BBSim');
    AAlT = AASim;
    ABlT = conj(ABSim); % * exp(-1i*3);
    BBlT = BBSim;
    %-------------------------------------------------------------------------%
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Real Data
    %-------------------------------------------------------------------------%
    load('AAlTermPZ');
    load('ABlTermPZ');
    load('BBlTermPZ');

%     ABlTerm = longTermAverage(AB,50);
%     AAlTerm = longTermAverage(AA,50);
%     BBlTerm = longTermAverage(BB,50);

    AAlT = AAlTerm;
    ABlT = ABlTerm; % .* exp(1i*peakPhaseInFrequencyCygA)'; % * conj(peakIntensityGain); %(sqrt(18) * 1024); % * exp(-1i*3);
    BBlT = BBlTerm; % BBlTerm / 18;
    %-------------------------------------------------------------------------%
end
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract Related Frequency and Time Range
%-------------------------------------------------------------------------%
AAlT = AAlT(fStart:fStop,tStart:end);
ABlT = ABlT(fStart:fStop,tStart:end);
BBlT = BBlT(fStart:fStop,tStart:end);
%-------------------------------------------------------------------------%

[fSize, tSize] = size(ABlT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna Coordinates
%-------------------------------------------------------------------------%
mainLatMap = 35.247227;
mainLongMap = -116.793272;

mainLat1Map = 35.247190;
mainLong1Map = -116.793287;
mainLat2Map = 35.247222;
mainLong2Map = -116.793316;
mainLat3Map = 35.247256;
mainLong3Map = -116.793300;
mainLat4Map = 35.247264;
mainLong4Map = -116.793251;
mainLat5Map = 35.247233;
mainLong5Map = -116.793221;
mainLat6Map = 35.247196;
mainLong6Map = -116.793239;

outLatMap = 35.247469;
outLongMap = -116.791509;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sky Source Coordinates
%-------------------------------------------------------------------------%
raCygnusA = 19 + 59/60 + 28.3/3600; % in Hours
decCygnusA = 40 + 44/60 + 02/3600; % in Degrees
% Cygnus A: 19h59m28.3s +40d44m02s
raCasA = 23 + 23/60 + 27.9/3600; % in Hours
decCasA = 58 + 48/60 + 42/3600; % in Degrees
% Cas A: 23h23m27.9s +58d48m42s
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the number of days from J2000.0 for 05:45:11 hrs UT on 2016 July 21th
%-------------------------------------------------------------------------%
startHour = 05 + 45/60 + 11/3600;
fractionDay = startHour / 24;
dayNum = 21;
daysSinceJanuary = 182;
daysSinceJ2000 = 5842.5;
numOfDays = daysSinceJ2000 + daysSinceJanuary + dayNum + fractionDay;

% Find the Local Siderial Time (LST) of Each Antenna
lstMain = 100.46 + 0.985647 * numOfDays + mainLongMap + 15*startHour;
lstMain = mod(lstMain,360);

lstMain1 = 100.46 + 0.985647 * numOfDays + mainLong1Map + 15*startHour;
lstMain1 = mod(lstMain1,360);
lstMain2 = 100.46 + 0.985647 * numOfDays + mainLong2Map + 15*startHour;
lstMain2 = mod(lstMain2,360);
lstMain3 = 100.46 + 0.985647 * numOfDays + mainLong3Map + 15*startHour;
lstMain3 = mod(lstMain3,360);
lstMain4 = 100.46 + 0.985647 * numOfDays + mainLong4Map + 15*startHour;
lstMain4 = mod(lstMain4,360);
lstMain5 = 100.46 + 0.985647 * numOfDays + mainLong5Map + 15*startHour;
lstMain5 = mod(lstMain5,360);
lstMain6 = 100.46 + 0.985647 * numOfDays + mainLong6Map + 15*startHour;
lstMain6 = mod(lstMain6,360);

lstOut = 100.46 + 0.985647 * numOfDays + outLongMap + 15*startHour;
lstOut = mod(lstOut,360);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
% nSAv = 50;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Total Scanned Hour Duration
%-------------------------------------------------------------------------%
h = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
h = h * pi/12; % in radians
hSize = length(h);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination (in degrees)
%-------------------------------------------------------------------------%
dec = linspace(0,80,50);
dSize = length(dec);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Right Ascention (in degrees and hours)
%-------------------------------------------------------------------------%
rAT = linspace(-90,0,100); % in degrees (global coordinates)
rA = (linspace(0,90,100) + lstMain)*12/180;
rSize = length(rA);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
fc = (0.09765625*fStart:0.09765625:0.09765625*fStop)'*1e6;
c = physconst('LightSpeed');
%-------------------------------------------------------------------------%


intensityBF = zeros(dSize,rSize);
intensityMF = zeros(dSize,rSize);
intensityCAPON = zeros(dSize,rSize);

detFactorR = 1 ./ (AAlT.*BBlT - abs(ABlT).^2);
detFactorR(isinf(detFactorR)) = 0;

% phaseFunctionCAPON = zeros(length(dec),length(ra),length(fc));
% phaseFunctionMF = zeros(length(dec),length(ra),length(fc));

if SimulationMode
    thao = 0;
else
    % nsF = 67; % 16; 68 ; 135; 69.5; % Best Value is 67
    % thao = nsF/(100e6);
    thao = 6.8677e-07;
end

for dD = 1:dSize
    for rR = 1:rSize
        
%         tgT = reshape(tgApproximate(dD,rR,:),1,length(h));
        tgT = reshape(tgExact(dD,rR,:,:),6,length(h)); % Time delay between outrigger and each LoFAMS ring antenna (6 by hour length)
        
        %%% Coherent
%         phaseComp = exp(-2*pi*1i*fc*(tgT-thao));
        phaseComp = zeros(length(fc),length(h));
        for nA = 1:6
            phaseComp = phaseComp + exp(-2*pi*1i*fc*(tgT(nA,:)-thao));
        end
        % phaseComp = phaseComp./abs(phaseComp);
        % phaseComp = phaseComp / 6;
        
        %%% Matched Filter
        intensityMF(dD,rR) = sum(sum(ABlT .* phaseComp));      

        %%% Beamforming
%         intensityBT = AAlT + 2*real(ABlT .* phaseComp) + BBlT;
%         intensityBF(dD,rR) = sum(sum(intensityBT));         
        
        %%% Capon
        intensityCT = detFactorR .* (AAlT - 2*real(ABlT .* phaseComp) + BBlT);
        
        % intensityCT = detFactorR .* (AAlT - 2*real(ABlT .* phaseComp) + BBlT.*abs(phaseComp).^2);
        
        % intensityCT = (2*abs(ABlT .* phaseComp)) ./ (AAlT+1e6);
        % intensityCT = detFactorR .* (- 2*real(ABlT .* phaseComp));
        % intensityCT = sum(intensityCT);
        intensityCT = 1./intensityCT;
        intensityCT(isinf(intensityCT)) = 0;
        intensityCAPON(dD,rR) = sum(sum(intensityCT));
%         
%         intensityMT = (AAlT + 2*real(ABlT.*phaseComp) + BBlT.*(abs(phaseComp).^2));
%         intensityMVDR(dD,rR) = sum(sum(intensityMT));
        
%         phaseFunctionMF(dD,rR,:) = sum(ABlT .* phaseComp , 2);
%         phaseFunctionCAPON(dD,rR,:) = sum(intensityCT , 2);
        
        %%% Nonceherent
%         phaseComp = exp(-2*pi*1i*fc*tgT);
%         phaseComp = 0;
%         for nA = 1:6
%             phaseComp = phaseComp + exp(-2*pi*1i*fc*tgT(nA,:));
%         end
%         
%         intensityMF(dD,rR) = sum(abs(sum(ABlT .* phaseComp , 2)).^2);
        
%         intensityBT = AAlT + 2*real(ABlT .* phaseComp) + BBlT;
%         intensityBF(dD,rR) = sum(abs(sum(intensityBT , 2)).^2);
%         
%         intensityCT = detFactorR .* (BBlT - 2*real(ABlT.*conj(phaseComp)) + AAlT.*(abs(phaseComp).^2));
%         intensityCAPON(dD,rR) = 1/sum(abs(sum(intensityCT , 2)).^2);
%         
%         intensityMT = (AAlT + 2*real(ABlT.*phaseComp) + BBlT.*(abs(phaseComp).^2));
%         intensityMVDR(dD,rR) = sum(abs(sum(intensityMT , 2)).^2);
        
    end
end

mesh(rA,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);

% figure
% mesh(rA,dec,abs(intensityBF),'FaceColor','interp','LineStyle','none')
% xlabel('RA (Hour)')
% ylabel('Dec (Degrees)')
% view(2)
% colormap('jet');
% set(gca,'xlim',[min(rA) max(rA)]);

figure
mesh(rA,dec,abs(intensityCAPON),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);

% figure
% mesh(ra,dec,abs(intensityMVDR))
% xlabel('RA (Degrees)')
% ylabel('Dec (Degrees)')

% CygAPhase = reshape(phaseFunctionMF(26,48,:),1,length(fc));
% figure
% plot(fc/1e6,unwrap(angle(CygAPhase)))
% xlabel('Frequency (MHz)')
% ylabel('Unwrapped Phase Angle (Degrees)')
% figure
% CasAPhase = reshape(phaseFunctionMF(38,77,:),1,length(fc));
% plot(fc/1e6,unwrap(angle(CasAPhase)))
% xlabel('Frequency (MHz)')
% ylabel('Unwrapped Phase Angle (Degrees)')