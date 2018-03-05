addpath(genpath('C:\Users\memin\Box Sync\LoFASMWork\Algorithms\CalibrationWork'))
addpath(genpath('C:\D\UTD\Research\LoFASM\ImagingWork'))
addpath(genpath('C:\D\UTD\Research\LoFASM\harddiskRecordedData'))

%-------------------------------------------------------------------------%
% Initial Variables
%-------------------------------------------------------------------------%
fRange = 250:750;
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

SimulationMode = true;
Whitening = false;
antennaGainFactor = false;

StationID = 4;

%-------------------------------------------------------------------------%
% Data Record Time
%-------------------------------------------------------------------------%
% Data Record Start Time [h,m,s,d,m,y]
dataRecordTime = [5,45,11,21,7,2016];

% dataRecordTime = [16,59,59,26,4,2017]; % start time
% dataRecordTime = [11,28,41,27,4,2017]; % ~18:00 LST
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
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
    ABlTerm = ABSim;
    %-------------------------------------------------------------------------%
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Real Data
    %-------------------------------------------------------------------------%
    if nSAv == 447
        load('ABlTermPZ');
        
        % load('AD_NS_18h_22h_lTerm');
        % ABlTerm = AD_NS_18h_22h_lTerm;
    else
        load('AB_Padded.mat')
        ABlTerm = longTermAverage(AB,nSAv);      
    end
end

if Whitening
    ABlTerm = whiteningApproach1(ABlTerm);
end

ABlT = ABlTerm; clear ABlTerm
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract Related Frequency and Time Range
%-------------------------------------------------------------------------%
ABlT = ABlT(fRange,:);
[fSize, tSize] = size(ABlT);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna Coordinates
%-------------------------------------------------------------------------%

if StationID == 1
    mainLat2Map = 26.555434;
    mainLong2Map = -97.442024;
    mainLat5Map = 26.555452;
    mainLong5Map = -97.441938;
    mainAltMap = 3;
    
    outLatMap = 26.556041;
    outLongMap = -97.440489;
    outAltMap = 3;
    
elseif StationID == 4
    % mainLat1Map = 35.247190;
    % mainLong1Map = -116.793287;
    mainLat2Map = 35.247222;
    mainLong2Map = -116.793316;
    % mainLat3Map = 35.247256;
    % mainLong3Map = -116.793300;
    % mainLat4Map = 35.247264;
    % mainLong4Map = -116.793251;
    mainLat5Map = 35.247233;
    mainLong5Map = -116.793221;
    % mainLat6Map = 35.247196;
    % mainLong6Map = -116.793239;
    mainAltMap = distdim(3547,'ft','m');
    
    outLatMap = 35.247469;
    outLongMap = -116.791509;
    uDistance = 3; % Zenith baseline in meters
    outAltMap = distdim(3525,'ft','m') + uDistance;
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sky Source Coordinates
%-------------------------------------------------------------------------%
raCygnusA = 19 + 59/60 + 28.3/3600; % in Hours (19.9912)
decCygnusA = 40 + 44/60 + 02/3600; % in Degrees (40.7339)
% Cygnus A: 19h59m28.3s +40d44m02s
raCasA = 23 + 23/60 + 27.9/3600; % in Hours (23.3911)
decCasA = 58 + 48/60 + 42/3600; % in Degrees (58.8117)
% Cas A: 23h23m27.9s +58d48m42s
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the LST
%-------------------------------------------------------------------------%
% LST = calculate_LST(longitude,hour,minute,second,day,month,year);
dataRecordTime = num2cell(dataRecordTime);
LST = calculate_LST(outLongMap,dataRecordTime{:});
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Total Scanned Hour Duration
%-------------------------------------------------------------------------%
h = (0 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baseline from Map Locations
%-------------------------------------------------------------------------%
wgs84 = wgs84Ellipsoid('meters');

arrayCenterLat = (mainLat2Map+mainLat5Map)/2;
arrayCenterLon = (mainLong2Map+mainLong5Map)/2;

[northMO,eastMO,downMO] = geodetic2ned(arrayCenterLat,arrayCenterLon,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);

innerRingRotation = 90-azimuth(mainLat2Map,mainLong2Map,mainLat5Map,mainLong5Map,wgs84); % Counter-clockwise - Degrees

rotNed2Body = [cosd(innerRingRotation) sind(innerRingRotation) 0; ...
                -sind(innerRingRotation) cosd(innerRingRotation) 0; ...
                0 0 1];
    
% rotNed2BodyMatlab = angle2dcm( (innerRingRotation)*pi/180, 0, 0 );

innerNED = [2.205*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
            2.205*[-1, -2, -1, 1, 2, 1];...
            zeros(1,6)];

% plot(innerNED(2,:),innerNED(1,:),'.')

innerNED = rotNed2Body * innerNED + [northMO;eastMO;downMO];

x = innerNED(2,:) + 0.3030; % dxEst;
y = innerNED(1,:); % + 0.6566; % dyEst;
z = -innerNED(3,:) + 2.1212; % dzEst;

% plot(x,y,'o') x: towards the east, y: toward the north
% sum((innerNED(:,1)-innerNED(:,2)).^2)^0.5
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Scanned Right Ascention and Total Scanned Hour Angle
%-------------------------------------------------------------------------%
rAStart = 18; rAStop = 24; % in hours
rA = rAStart: (nSAv*0.083886)/3600: rAStop;

H = LST - rA; H = flip(H);
H = [H H(end)+h(2:end)];
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination
%-------------------------------------------------------------------------%
dec = linspace(30,70,300).'; % in degrees
% dec = decCygnusA;
% dec = decCasA;
dSize = length(dec);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baselines and Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = x; X = reshape(X,1,1,6);
Y = z*sind(l) + y*cosd(l); Y = reshape(Y,1,1,6);
Z = z*cosd(l) - y*sind(l); Z = reshape(Z,1,1,6);

w = -cosd(dec)*sin(H*pi/12).*X + sind(dec).*Y + cosd(dec)*cos(H*pi/12).*Z;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
deltaF = 0.09765625*1e6; % 100e6/1024
indF = 0:1023;
fc = indF*deltaF;
fc = fc(fRange);
c = physconst('LightSpeed');
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Instrument Delay
%-------------------------------------------------------------------------%
if SimulationMode
    tI = -6.9000e-07;
else
    % tI = 67/(100e6);
    tI = -6.9000e-07; % tIEst;
    % tI = -6.8728e-07; % 6.6055e-07; % 6.74609375e-07; % 6.8677e-07; % 68.73/(100e6); % ones(1,6)*6.7000e-07; % 4.0371e-07; %  % 0.0068; % -0.0011; %
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Exact Matched Filter with 6 Main Antenna is Used
%---------------------------------------------------------------------%
if fSize>1
    fc = reshape(fc,1,1,fSize);
end

[~,antGain] = radiationPatternLofasm(fc,H);
antGain = repmat(antGain,[dSize 1 1]);
antGain = cosd(dec-outLatMap) .* antGain;

matchedFilter = zeros(dSize,length(H),fSize);
matchedFilterEachAntenna = zeros(dSize,length(H),fSize);
for nA = 1:6
    tG = squeeze(w(:,:,nA))/c;  
    if length(tI) == 6
        matchedFilterEachAntenna = exp(-1i*2*pi*fc.*(tG+tI(nA)));
    elseif length(tI) == 1
        matchedFilterEachAntenna = exp(-1i*2*pi*fc.*(tG+tI)); 
    end
    matchedFilter = matchedFilter + matchedFilterEachAntenna;
end

if antennaGainFactor
    matchedFilter = antGain .* matchedFilter;
end
%---------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare AB for IFFT
%-------------------------------------------------------------------------%
ABlT = [ABlT zeros(fSize,length(H)-tSize)];
if fSize>1
    ABlT = reshape(transpose(ABlT),1,length(H),fSize);
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D FFT
if fSize==1
    ABlTfft = fft(ABlT);
else
    ABlTfft = fft(ABlT,length(H),2);
end
matchedFilterfft = fft2(flip(matchedFilter,2));
%-------------------------------------------------------------------------%

intensityMFT = sum(flip(ifft2(ABlTfft.*matchedFilterfft),2),3);
intensityMF = flip(intensityMFT(:,1:length(rA)),2);

if (dSize == 1)
    plot(rA,abs(intensityMF))
    xlabel('RA (Hour)')
else
    mesh(rA,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
    xlabel('RA (Hour)')
    ylabel('Dec (Degrees)')
    view(2)
    colormap('jet');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the LST
%-------------------------------------------------------------------------%
function LST = calculate_LST(longitude,hour,minute,second,day,month,year)

% JD2000 = juliandate(2000,1,1,12,0,0);
JD2000 = juliandate(2000,1,1,11,58,55.816);

JD = juliandate(year,month,day);
UT = hour + minute/60 + second/3600;
T = (JD - JD2000)/36525.0;
GST = 6.697374558 + (2400.051336*T)+(0.000025862*T^2)+(UT*1.0027379093);
GST = mod(GST,24);

% Find the Local Siderial Time (LST)
LST = mod(longitude/15 + GST,24);
end
%-------------------------------------------------------------------------%

%%% Test FFT
% x = [1 2 8 16];
% y = [1 2 6 8 1 2 5 7];
% z1 = x(1)*y(1:5) + x(2)*y(2:6)+x(3)*y(3:7)+x(4)*y(4:8);
% x2 = [x zeros(1,length(y)-length(x))];
% z2 = flip(ifft(fft(x2).*fft(flip(y))));
% z2 = z2(1:5);

% x = ABlT;
% y = matchedFilter;
% 
% x2 = [x zeros(1,length(y)-length(x))];
% z2 = flip(ifft(fft(x2,length(y)).*fft(flip(y))));
% z2 = flip(z2(1:length(rA)));
% 
% z = zeros(1,length(rA));
% for n = 1:tSize
%    z = z + x(n)*y(n:length(rA)+n-1);
% end