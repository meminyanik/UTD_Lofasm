addpath(genpath('C:\D\UTD\Research\LoFASM\CalibrationWork'))
addpath(genpath('C:\D\UTD\Research\LoFASM\harddiskRecordedData'));
addpath(genpath('C:\D\UTD\Research\LoFASM\ImagingWork'));

%-------------------------------------------------------------------------%
% Initial Variables
%-------------------------------------------------------------------------%
fRange = 549:550;
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

SimulationMode = true;
Whitening = false;

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
% nSAv = 1;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Real Data
%-------------------------------------------------------------------------%
if (~SimulationMode)
    if nSAv == 447
        load('AAlTermPZ');
        load('ABlTermPZ');
        load('BBlTermPZ');
    else
        load('AA_Padded.mat')
        load('AB_Padded.mat')
        load('BB_Padded.mat')
        AAlTerm = longTermAverage(AA,nSAv);
        ABlTerm = longTermAverage(AB,nSAv);
        BBlTerm = longTermAverage(BB,nSAv);
    end
    
    if Whitening
        ABlTerm = whiteningApproach1(ABlTerm);
    end
    
    AAlT = AAlTerm; clear AAlTerm
    ABlT = ABlTerm; clear ABlTerm
    BBlT = BBlTerm; clear BBlTerm

else
    AAlT = AASim;
    ABlT = ABSim;
    BBlT = BBSim;
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract Related Frequency and Time Range
%-------------------------------------------------------------------------%
AAlT = AAlT(fRange,:);
ABlT = ABlT(fRange,:);
BBlT = BBlT(fRange,:);
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

% plot(innerNED(2,:),innerNED(1,:),'o')

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
rA = linspace(18,24,200); % in hours
% rA = raCygnusA;
% rA = raCasA;
rSize = length(rA);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination
%-------------------------------------------------------------------------%
dec = linspace(30,70,100).'; % in degrees
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
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
deltaF = 0.09765625*1e6; % 100e6/1024
indF = 0:1023;
fc = indF*deltaF;
fc = fc(fRange).';
c = physconst('LightSpeed');
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate Determinant Factor
%-------------------------------------------------------------------------%
detFactorR = 1 ./ (AAlT.*BBlT+1e-5 - abs(ABlT).^2);
detFactorR(isinf(detFactorR)) = 0;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Instrument Delay
%-------------------------------------------------------------------------%
% tI = 67/(100e6);
tI = -6.9000e-07; % tIEst;
tI = repmat(tI,1,1,6);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize Sky Intensity
%-------------------------------------------------------------------------%
intensityBF = zeros(dSize,rSize);
intensityMF = zeros(dSize,rSize);
intensityCAPON = zeros(dSize,rSize);
%-------------------------------------------------------------------------%


for dD = 1:dSize
    for rR = 1:rSize
        
        H = (LST - rA(rR) + h)*pi/12; % Hour Angle in radians
        w = -cosd(dec(dD))*sin(H).*X + sind(dec(dD)).*Y + cosd(dec(dD))*cos(H).*Z;
        tG = w/c;
        matchedFilter = sum(exp(-1i*2*pi*fc.*(tG + tI)),3);
        
        %%% Matched Filter
        intensityMF(dD,rR) = sum(sum(ABlT .* matchedFilter));      

        %%% Beamforming
        intensityBT = AAlT + 2*real(ABlT .* matchedFilter) + BBlT;
        intensityBF(dD,rR) = sum(sum(intensityBT));         
        
        %%% Capon
        % intensityCT = detFactorR .* (AAlT - 2*real(ABlT .* matchedFilter) + BBlT);
        intensityCT = AAlT - 2*real(ABlT .* matchedFilter) + BBlT;
        % intensityCT = (AAlT - 2*real(ABlT .* matchedFilter) + BBlT.*abs(matchedFilter).^2);
        % intensityCT = (2*abs(ABlT .* phaseComp)) ./ (AAlT+1e6);
        % intensityCT = detFactorR .* (- 2*real(ABlT .* phaseComp));
        % intensityCT = sum(intensityCT);
        intensityCT = 1./intensityCT;
        intensityCT(isinf(intensityCT)) = 0;
        intensityCAPON(dD,rR) = sum(sum(intensityCT));
        
    end
end

mesh(rA,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Processing on Beamforming
%-------------------------------------------------------------------------%
intensityBF_P = intensityBF - mean(mean(intensityBF));
intensityBF_P = (hilbert(intensityBF_P.')).';

figure
mesh(rA,dec,abs(intensityBF_P),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post Processing on CAPON
%-------------------------------------------------------------------------%
intensityCAPON_P = intensityCAPON - mean(mean(intensityCAPON));
intensityCAPON_P = (hilbert(intensityCAPON_P.')).';

figure
mesh(rA,dec,abs(intensityCAPON_P),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);

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