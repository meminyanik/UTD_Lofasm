tStart = 1;
fStart = 110; % 140; % 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz
fStop = 550; %800; % 1024 = 100 MHz, 820 = 80 MHz, 417 = 40 MHz

SimulationMode = false;

load('tgExactConvolvedH_LST_V2.mat')
% load('ABFilteredlTerm_nSample120_Dist1e34.mat')

load('gainFunctionConvolvedH_LST_V2.mat')

load('peakPhaseInFrequencyCygAwhitening');

% peakPhaseInFrequencyCygA = unwrap(peakPhaseInFrequencyCygA);
% load('peakPhaseInTimeCygA');

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
    % load('ABSim-CygA');
    ABlT = ABSim;
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
% Find the LST for 05:45:11 hrs UT on 2016 July 21th
%-------------------------------------------------------------------------%
hour = 5; minute = 45; second = 11;
day = 21; month = 7; year = 2016;

% JD2000 = juliandate(2000,1,1,12,0,0);
JD2000 = juliandate(2000,1,1,11,58,55.816);

JD = juliandate(year,month,day);
UT = hour + minute/60 + second/3600;
T = (JD - JD2000)/36525.0;
GST = 6.697374558 + (2400.051336*T)+(0.000025862*T^2)+(UT*1.0027379093);
GST = mod(GST,24);

% Find the Local Siderial Time (LST) of Each Antenna

lstMain1 = mod(mainLong1Map + GST*15,360);
lstMain2 = mod(mainLong2Map + GST*15,360);
lstMain3 = mod(mainLong3Map + GST*15,360);
lstMain4 = mod(mainLong4Map + GST*15,360);
lstMain5 = mod(mainLong5Map + GST*15,360);
lstMain6 = mod(mainLong6Map + GST*15,360);

lstOut = mod(outLongMap + GST*15,360);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Total Scanned Hour Duration
%-------------------------------------------------------------------------%
raLd = 0; raRd = 90;

hT = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
raLs = round(raLd*12/180*3600/(nSAv*0.083886));
raRs = round(raRd*12/180*3600/(nSAv*0.083886));

raHl = (raLs*nSAv*0.083886 : nSAv*0.083886 : -1*nSAv*0.083886)/3600;
raHr = (tSize*nSAv*0.083886 : nSAv*0.083886 : (tSize+raRs)*nSAv*0.083886)/3600;
h = [raHl hT raHr] * pi/12; % in radians
tSizeN = length(h);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination (in degrees)
%-------------------------------------------------------------------------%
dec = linspace(0,80,50);
dSize = length(dec);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Right Ascention (in hours)
%-------------------------------------------------------------------------%
rA = (raLs*nSAv*0.083886 : nSAv*0.083886 : raRs*nSAv*0.083886)/3600;
rA = rA + lstOut*12/180;
rSize = length(rA);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
fc = (0.09765625*(fStart-1):0.09765625:0.09765625*(fStop-1))'*1e6;
c = physconst('LightSpeed');
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare AB for IFFT
%-------------------------------------------------------------------------%
[~, padSl] = size(raHl);
[~, padSr] = size(raHr);
ABlT = [zeros(fSize,padSl) ABlT zeros(fSize,padSr)];
ABlT = fliplr(ABlT);
%-------------------------------------------------------------------------%

ABlT = reshape(transpose(ABlT),1,tSizeN,fSize);
fc = reshape(fc,1,1,fSize);
% If phase correction is done at Matched Filter
% peakPhaseInFrequencyCygA = reshape(peakPhaseInFrequencyCygA,1,1,fSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D FFT
ABlTfft = fft2(ABlT,dSize,tSizeN);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Instrument Delay
%-------------------------------------------------------------------------%
if SimulationMode
    tI = 0;
else
%     nsF = 67; % 16; 68 ; 135; 69.5; % Best Value is 67
%     tI = nsF/(100e6)*ones(1,6);
    tI = ones(1,6)*7.1845e-06;
    % tI = 0.0068 * ones(1,6); % T_newton; % ones(1,6)*6.7000e-07; % 4.0371e-07; %  % 0.0068; % -0.0011; % 
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Approximate Matched Filter with 1 Main Antenna is Used
%---------------------------------------------------------------------%
% phaseComp = exp(-2*pi*1i*fc.*(tgApproximate-thao));
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Exact Matched Filter with 6 Main Antenna is Used
%---------------------------------------------------------------------%
phaseComp = 0;
for nA = 1:6
    tgT = reshape(tgExact(:,nA,:),dSize,tSizeN);
    phaseComp = phaseComp + exp(-2*pi*1i*fc.*(tgT-tI(nA))); % .* exp(-1i*peakPhaseInFrequencyCygA); % If phase correction is done at Data;
end
% gainFunction = repmat(gainFunction,[1,1,fSize]);
% phaseComp = phaseComp .* gainFunction;
%---------------------------------------------------------------------%

intensityMFT = sum(ifft2(ABlTfft.*fft2(phaseComp)),3);
intensityMF = intensityMFT(:,1:padSr);
intensityMF = fliplr(intensityMF);
% Old Version
% intensityMF = circshift(intensityMFT,padSl,2);
% intensityMF = fliplr(intensityMF);
% intensityMF = intensityMF(:,length(hT)+1:end);

% mesh(rA,dec,pow2db(abs(intensityMF)),'FaceColor','interp','LineStyle','none')
mesh(rA,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);

% contour(rA,dec,abs(intensityMF))
% xlabel('RA (Hour)')
% ylabel('Dec (Degrees)')
% view(2)
% colormap('gray');
% set(gca,'xlim',[min(rA) max(rA)]);

% figure
% mesh(abs(intensityMFT),'FaceColor','interp','LineStyle','none')
% view(2)
% colormap('jet');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Peak Phase Analysis
%---------------------------------------------------------------------%
% intensityMFT = ifft2(ABlTfft.*fft2(phaseComp));
% 
% intensity = zeros(50,577);
% 
% peakPhaseInFrequencyCygA = zeros(1,length(fc));
% peakPhaseInFrequencyCasA = zeros(1,length(fc));
% 
% figure
% for nF = 1:length(fc)
%     intensityMF = intensityMFT(:,:,nF);
%     intensityMF = circshift(intensityMF,padSl,2);
%     intensityMF = fliplr(intensityMF);
%     intensityMF = intensityMF(:,length(hT)+1:end);
%     
%     peakPhaseInFrequencyCygA(nF) = angle(intensityMF(26,182));
%     peakPhaseInFrequencyCasA(nF) = angle(intensityMF(37,514));
%     
% %     intensity = intensity + intensityMF;
%     
% 
% %     mesh(rA,dec,abs(intensity),'FaceColor','interp','LineStyle','none')
% %     xlabel('RA (Hour)')
% %     ylabel('Dec (Degrees)')
% %     view(2)
% %     colormap('jet');
% %     set(gca,'xlim',[min(rA) max(rA)]);
% %     drawnow
%     
% % %     plot(unwrap(peakPhase));
% %     plot(peakPhaseInFrequency);
% %     drawnow;
%     
% end
% 
% nsF = 15; % 68 ; 135; 69.5;
% thao = nsF/(100e6);
% fc = (0.09765625*fStart:0.09765625:0.09765625*fStop)'*1e6;
% 
% figure
% plot(angle(exp(1i*peakPhaseInFrequencyCygA)))
% hold on
% plot(angle(exp(-2*pi*1i*fc*thao)))
% 
% figure
% plot(unwrap(peakPhaseInFrequencyCygA))
% hold on
% plot(unwrap(-2*pi*fc*thao))
%---------------------------------------------------------------------%