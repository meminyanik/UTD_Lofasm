%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
currentPath = pwd();
addpath(genpath([currentPath,'/../../../Algorithms']));
addpath(genpath([currentPath,'/../../../RecordedData']));

%-------------------------------------------------------------------------%
% Initial Variables
%-------------------------------------------------------------------------%
fRange = 550:550;
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

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
% nSAv = 447;
% nSAv = 50;
nSAv = 1;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%-------------------------------------------------------------------------%
if nSAv == 447
    % load('ABlTermPZ');
    load('AD_NS_18h_22h_lTerm');
    ABlTerm = AD_NS_18h_22h_lTerm;
elseif nSAv == 1
    % load('AB_Padded.mat')
    ABlTerm = AB;
else
    load('AB_Padded.mat')
    ABlTerm = longTermAverage(AB,nSAv);
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

x = innerNED(2,:);
y = innerNED(1,:);
z = @(dz) -innerNED(3,:) + dz; %- 23.0; % 4.5; %% repmat(3,1,6); % % repmat(0.6907,1,6); %

% plot(x,y,'o') x: towards the east, y: toward the north
% sum((innerNED(:,1)-innerNED(:,2)).^2)^0.5
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baselines and Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = x;
Y = @(z) z*sind(l) + y*cosd(l);
Z = @(z) z*cosd(l) - y*sind(l);

dCygA = decCygnusA;
rCygA = raCygnusA;
HCygA = (LST + h - rCygA).';
wCygA = @(z) -cosd(dCygA)*sin(HCygA*pi/12)*X + sind(dCygA)*Y(z) + cosd(dCygA)*cos(HCygA*pi/12)*Z(z);

dCasA = decCasA;
rCasA = raCasA;
HCasA = (LST + h - rCasA).';
wCasA = @(z) -cosd(dCasA)*sin(HCasA*pi/12)*X + sind(dCasA)*Y(z) + cosd(dCasA)*cos(HCasA*pi/12)*Z(z);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Exact Matched Filter with 6 Main Antenna is Used
%---------------------------------------------------------------------%
% gainCygA = repmat(cosd(HCygA*15).',fSize,1);
[~,gainCygA] = radiationPatternLofasm(fc,HCygA.');
tGcygA = @(z) reshape(wCygA(z),1,tSize,6)/c;
matchedFilterCygA = @(z) sum(exp(-1i*2*pi*fc.*tGcygA(z)),3);
if antennaGainFactor
    matchedFilterCygA = @(z) gainCygA .* matchedFilterCygA(z);
end

% gainCasA = repmat(cosd(HCasA*15).',fSize,1);
[~,gainCasA] = radiationPatternLofasm(fc,HCasA.');
tGcasA = @(z) reshape(wCasA(z),1,tSize,6)/c;
matchedFilterCasA = @(z) sum(exp(-1i*2*pi*fc.*tGcasA(z)),3);
if antennaGainFactor
    matchedFilterCasA = @(z) gainCasA .* matchedFilterCasA(z);
end

% Simulate AB Data
% ABSim = conj(matchedFilterCygA)+conj(matchedFilterCasA);
%---------------------------------------------------------------------%


%---------------------------------------------------------------------%
cygAPeakIntensity = @(z) sum(ABlT.*matchedFilterCygA(z),2);
% figure; plot(fc,unwrap(angle(cygAPeakIntensity)));
casAPeakIntensity = @(z) sum(ABlT.*matchedFilterCasA(z),2);
% figure; plot(fc,unwrap(angle(casAPeakIntensity)));


nFFT = 1*1024;
cygAPeakIntensityFFT = @(z) fftshift(fft(cygAPeakIntensity(z),nFFT));
casAPeakIntensityFFT = @(z) fftshift(fft(casAPeakIntensity(z),nFFT));

estimatedZrange = linspace(-30,20,100);
cygAPeakIntensityFFTallZ = zeros(nFFT,length(estimatedZrange));
casAPeakIntensityFFTallZ = zeros(nFFT,length(estimatedZrange));
for nZ = 1:length(estimatedZrange)
    cygAPeakIntensityFFTallZ(:,nZ) = cygAPeakIntensityFFT(estimatedZrange(nZ));
    casAPeakIntensityFFTallZ(:,nZ) = casAPeakIntensityFFT(estimatedZrange(nZ));
end
mesh(estimatedZrange,0:1023,abs(cygAPeakIntensityFFTallZ))
figure
mesh(estimatedZrange,0:1023,abs(casAPeakIntensityFFTallZ))

[~,indTiCygA] = max(abs(cygAPeakIntensityFFTallZ));
[~,indTiCasA] = max(abs(casAPeakIntensityFFTallZ));

tIEstCygA = (indTiCygA-1-nFFT/2)/(nFFT*deltaF);
tIEstCasA = (indTiCasA-1-nFFT/2)/(nFFT*deltaF);

figure; plot(indTiCygA)
hold on
plot(indTiCasA)


% figure; plot(abs(cygAPeakIntensityFFT));
% figure; plot(abs(casAPeakIntensityFFT));


%---------------------------------------------------------------------%


%---------------------------------------------------------------------%

% deltaZ = c*(tIEstCygA - tIEstCasA)/(sind(l)*(sind(dCygA)-sind(dCasA)));

%---------------------------------------------------------------------%


%---------------------------------------------------------------------%
% cygAPeakIntensity2D = ABlT.*matchedFilterCygA;
% casAPeakIntensity2D = ABlT.*matchedFilterCasA;
% 
% cygAPeakIntensity2DFFT = ifft(cygAPeakIntensity2D,nFFT/128);
% casAPeakIntensity2DFFT = ifft(casAPeakIntensity2D,nFFT/128);
% 
% [~,indTiCygA2D] = max(abs(cygAPeakIntensity2DFFT));
% [~,indTiCasA2D] = max(abs(casAPeakIntensity2DFFT));
% 
% tIEstCygA2D = (indTiCygA2D-1-nFFT/2)/(nFFT*deltaF);
% tIEstCasA2D = (indTiCasA2D-1-nFFT/2)/(nFFT*deltaF);

% tIAnalyticCygA = @(tI,dZ) tI + dZ*(sind(dCygA)*sind(l)+cosd(dCygA)*cosd(l)*cos(HCygA.'*pi/12));
% tIAnalyticCasA = @(tI,dZ) tI + dZ*(sind(dCasA)*sind(l)+cosd(dCasA)*cosd(l)*cos(HCasA.'*pi/12));
% 
% tIEstCygA = (444-1-nFFT/2)/(nFFT*deltaF);
% modelCygA = exp(1i*2*pi*indF.'*deltaF*tIAnalyticCygA(tIEstCygA,10));
%---------------------------------------------------------------------%


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


