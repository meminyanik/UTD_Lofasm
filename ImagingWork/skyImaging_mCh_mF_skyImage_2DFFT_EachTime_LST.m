tStart = 1;
fStart = 550; % 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz
fStop = 550; % 1024 = 100 MHz, 820 = 80 MHz, 417 = 40 MHz

SimulationMode = false;

% load('tgExact_LST_V2.mat')

% load('peakPhaseInFrequencyCasAwhitening');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
% nSAv = 50;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%-------------------------------------------------------------------------%
if SimulationMode
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Simulated Data
    %-------------------------------------------------------------------------%
    load('ABSimCygA');
    ABlT = ABSimCygA;
    % ABlT = conj(ABSimCygA);
    %-------------------------------------------------------------------------%
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Real Data
    %-------------------------------------------------------------------------%
    if nSAv == 447
        % load('AAlTermPM');
        % AAlT = AAlTerm;
        
        load('ABlTermPZ');
        ABlTerm = whiteningApproach1(ABlTerm);
        ABlT = ABlTerm; % .* exp(1i*peakPhaseInFrequencyCygA)'; % If phase correction is done at Data;
        % ABlT = ABlTerm;
        
        % load('BBlTermPM');
        % BBlT = BBlTerm; 
        
        % ABlT = ABFiltered;
 
    else
        % load('AA_Padded.mat')
        % AAlTerm = longTermAverage(AA,nSAv);
        % AAlT = AAlTerm;
        
        % load('AB_Padded.mat')
        % ABlTerm = longTermAverage(AB,nSAv);
        % ABlT = ABlTerm / 1024;
        
        % load('BB_Padded.mat')
        % BBlTerm = longTermAverage(BB,nSAv);
        % BBlT = BBlTerm;
    end
end
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract Related Frequency and Time Range
%-------------------------------------------------------------------------%
% AAlT = AAlT(fStart:fStop,tStart:end);
ABlT = ABlT(fStart:fStop,tStart:end);
% BBlT = BBlT(fStart:fStop,tStart:end);
[fSize, tSize] = size(ABlT);
%-------------------------------------------------------------------------%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Antenna Coordinates
% %-------------------------------------------------------------------------%
% mainLatMap = 35.247227;
% mainLongMap = -116.793272;
% 
% mainLat1Map = 35.247190;
% mainLong1Map = -116.793287;
% mainLat2Map = 35.247222;
% mainLong2Map = -116.793316;
% mainLat3Map = 35.247256;
% mainLong3Map = -116.793300;
% mainLat4Map = 35.247264;
% mainLong4Map = -116.793251;
% mainLat5Map = 35.247233;
% mainLong5Map = -116.793221;
% mainLat6Map = 35.247196;
% mainLong6Map = -116.793239;
% 
% outLatMap = 35.247469;
% outLongMap = -116.791509;
% %-------------------------------------------------------------------------%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Sky Source Coordinates
% %-------------------------------------------------------------------------%
% raCygnusA = 19 + 59/60 + 28.3/3600; % in Hours
% decCygnusA = 40 + 44/60 + 02/3600; % in Degrees
% % Cygnus A: 19h59m28.3s +40d44m02s
% raCasA = 23 + 23/60 + 27.9/3600; % in Hours
% decCasA = 58 + 48/60 + 42/3600; % in Degrees
% % Cas A: 23h23m27.9s +58d48m42s
% %-------------------------------------------------------------------------%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Find the number of days from J2000.0 for 05:45:11 hrs UT on 2016 July 21th
% %-------------------------------------------------------------------------%
% startHour = 05 + 45/60 + 11/3600;
% fractionDay = startHour / 24;
% dayNum = 21;
% daysSinceJanuary = 182;
% daysSinceJ2000 = 5842.5;
% numOfDays = daysSinceJ2000 + daysSinceJanuary + dayNum + fractionDay;
% 
% % Find the Local Siderial Time (LST) of Each Antenna
% lstMain = 100.46 + 0.985647 * numOfDays + mainLongMap + 15*startHour;
% lstMain = mod(lstMain,360);
% 
% lstMain1 = 100.46 + 0.985647 * numOfDays + mainLong1Map + 15*startHour;
% lstMain1 = mod(lstMain1,360);
% lstMain2 = 100.46 + 0.985647 * numOfDays + mainLong2Map + 15*startHour;
% lstMain2 = mod(lstMain2,360);
% lstMain3 = 100.46 + 0.985647 * numOfDays + mainLong3Map + 15*startHour;
% lstMain3 = mod(lstMain3,360);
% lstMain4 = 100.46 + 0.985647 * numOfDays + mainLong4Map + 15*startHour;
% lstMain4 = mod(lstMain4,360);
% lstMain5 = 100.46 + 0.985647 * numOfDays + mainLong5Map + 15*startHour;
% lstMain5 = mod(lstMain5,360);
% lstMain6 = 100.46 + 0.985647 * numOfDays + mainLong6Map + 15*startHour;
% lstMain6 = mod(lstMain6,360);
% 
% lstOut = 100.46 + 0.985647 * numOfDays + outLongMap + 15*startHour;
% lstOut = mod(lstOut,360);
% %-------------------------------------------------------------------------%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Time Axis
% %-------------------------------------------------------------------------%
% h = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
% h = (h + lstMain*12/180) * pi/12;
% %-------------------------------------------------------------------------%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Calculate Total Scanned Hour Duration
% %-------------------------------------------------------------------------%
% raLd = 0; raRd = 90;
% 
% hT = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
% raLs = round(raLd*12/180*3600/(nSAv*0.083886));
% raRs = round(raRd*12/180*3600/(nSAv*0.083886));
% 
% raHl = (raLs*nSAv*0.083886 : nSAv*0.083886 : -1*nSAv*0.083886)/3600;
% raHr = (tSize*nSAv*0.083886 : nSAv*0.083886 : (tSize+raRs)*nSAv*0.083886)/3600;
% h = [raHl hT raHr] * pi/12; % in radians
% tSizeN = length(h);
% %-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination (in degrees)
%-------------------------------------------------------------------------%
DEC = linspace(30,70,200);
dSize = length(DEC);
%-------------------------------------------------------------------------%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Define Right Ascention (in hours)
% %-------------------------------------------------------------------------%
% rA = (raLs*nSAv*0.083886 : nSAv*0.083886 : raRs*nSAv*0.083886)/3600;
% rA = rA + lstMain*12/180;
% rSize = length(rA);
% %-------------------------------------------------------------------------%
RA = linspace(18,24,300); % in hours
rSize = length(RA);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
fc = (0.09765625*(fStart-1):0.09765625:0.09765625*(fStop-1))'*1e6;
c = physconst('LightSpeed');
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare AB for IFFT
%-------------------------------------------------------------------------%
ABlT = reshape(transpose(ABlT),1,tSize,fSize);
fc = reshape(fc,1,1,fSize);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare Intensity
%-------------------------------------------------------------------------%
intensityMF = zeros(dSize,rSize);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Instrument Delay
%-------------------------------------------------------------------------%
if SimulationMode
    thao = 0;
else
    nsF = 67; % 16; 68 ; 135; 69.5; % Best Value is 67
    thao = nsF/(100e6);
%     thao = 0;
end
%-------------------------------------------------------------------------%

peakPhaseInTime = zeros(1,tSize);

figure
xlabel('RA (Degrees)')
ylabel('Dec (Degrees)')
colormap('jet');

for hH = 1:tSize
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% When Exact Matched Filter with 6 Main Antenna is Used
    %---------------------------------------------------------------------%
    tgT = reshape(tgExact(:,:,:,hH),dSize,rSize,6);
    phaseComp = 0;
    for nA = 1:6
        tgTA = reshape(tgT(:,:,nA),dSize,rSize);
        phaseComp = phaseComp + exp(-2*pi*1i*fc.*(tgTA-thao));
    end
    %---------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 2D FFT of AB Data
    ABlTfft = fft2(ABlT(:,hH,:),dSize,rSize);
    %-------------------------------------------------------------------------%
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate Intensity
    %-------------------------------------------------------------------------%
    intensityMF = intensityMF + sum(ifft2(ABlTfft.*fft2(phaseComp)),3);
    %-------------------------------------------------------------------------%
    
    mesh(RA,DEC,abs(intensityMF),'FaceColor','interp','LineStyle','none')
    view(2)
    drawnow;
    
%     peakPhaseInTime(hH) = angle(intensityMF(26,29));
%     plot(unwrap(peakPhaseInTime));
%     drawnow;

end

% plot(unwrap(peakPhaseInTime));

% figure
% mesh(ra,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
% xlabel('RA (Degrees)')
% ylabel('Dec (Degrees)')
% colormap('jet');
% view(2)