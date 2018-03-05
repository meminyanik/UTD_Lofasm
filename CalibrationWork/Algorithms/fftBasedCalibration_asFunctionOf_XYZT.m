addpath(genpath('C:\D\UTD\Research\LoFASM\CalibrationWork'))

%% Initial Variables
%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
currentPath = pwd();
addpath(genpath([currentPath,'/../../../Algorithms']));
addpath(genpath([currentPath,'/../../../RecordedData']));

fRange = 300:700;
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

SimulationMode = false;
Whitening = true;
antennaGainFactor = false;

SourceModelMode = 1; % 1: CygA, 2: Cas A, 3: CygA+CasA

StationID = 1;
arrayConfig = 0;
% 0: inner ring - outtrigger
% 1: outer ring - outtrigger
% 2: inner ring - outer ring
%-------------------------------------------------------------------------%

%% Data Record Time
% Data Record Start Time [h,m,s,d,m,y]
if (StationID == 1)
    dataRecordTime = [20,30,20,3,11,2017]; % 16.89 LST
elseif (StationID == 4)
    dataRecordTime = [5,45,11,21,7,2016];
end
%-------------------------------------------------------------------------%

%% Number of Samples Used During Averaging
% nSAv = 447; % Averaged value
nSAv = 200; % Whitened then averaged value
%-------------------------------------------------------------------------%

%% Read Correlation Data
readCorrelationData;
%-------------------------------------------------------------------------%

%% Subtract Related Frequency and Time Range
RlT = RlT(fRange,:);
[fSize, tSize] = size(RlT);
%-------------------------------------------------------------------------%

%% Get Antenna Coordinates
getAntennaCoordinates;
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

arrayCenterLat = (inMainLat2Map+inMainLat5Map)/2;
arrayCenterLon = (inMainLong2Map+inMainLong5Map)/2;

[northMO,eastMO,downMO] = geodetic2ned(arrayCenterLat,arrayCenterLon,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);

innerRingRotation = 90-azimuth(inMainLat2Map,inMainLong2Map,inMainLat5Map,inMainLong5Map,wgs84); % Counter-clockwise - Degrees
% outerRingRotation = 90-azimuth(outMainLat2Map,outMainLong2Map,outMainLat5Map,outMainLong5Map,wgs84); % Counter-clockwise - Degrees
outerRingRotation = innerRingRotation + 30;

rotInNed2Body = [cosd(innerRingRotation) sind(innerRingRotation) 0; ...
    -sind(innerRingRotation) cosd(innerRingRotation) 0; ...
    0 0 1];

rotOutNed2Body = [cosd(outerRingRotation) sind(outerRingRotation) 0; ...
    -sind(outerRingRotation) cosd(outerRingRotation) 0; ...
    0 0 1];

% rotNed2BodyMatlab = angle2dcm( (innerRingRotation)*pi/180, 0, 0 );

innerNED = [4.41/2*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
    4.41/2*[-1, -2, -1, 1, 2, 1];...
    zeros(1,6)];

outerNED = [7.6383/2*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
    7.6383/2*[-1, -2, -1, 1, 2, 1];...
    zeros(1,6)];

% plot(innerNED(2,:),innerNED(1,:),'.')

innerNED = rotInNed2Body * innerNED + [northMO;eastMO;downMO];
outerNED = rotOutNed2Body * outerNED + [northMO;eastMO;downMO];

if (arrayConfig == 0)
    x = @(dx) innerNED(2,:) + dx;
    y = @(dy) innerNED(1,:) + dy;
    z = @(dz) -innerNED(3,:) + dz;
elseif (arrayConfig == 1)
    x = @(dx) outerNED(2,:) + dx;
    y = @(dy) outerNED(1,:) + dy;
    z = @(dz) -outerNED(3,:) + dz;
end

plot(x(0),y(0),'o')
% hold on
% plot(x(0),y_i_ot(0),'o')
% x: towards the east, y: toward the north
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baselines and Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = @(dx) x(dx);
Y = @(dy,dz) z(dz)*sind(l) + y(dy)*cosd(l);
Z = @(dy,dz) z(dz)*cosd(l) - y(dy)*sind(l);

dCygA = decCygnusA;
rCygA = raCygnusA;
HCygA = (LST + h - rCygA).';
wCygA = @(dx,dy,dz) -cosd(dCygA)*sin(HCygA*pi/12)*X(dx) + sind(dCygA)*Y(dy,dz) + cosd(dCygA)*cos(HCygA*pi/12)*Z(dy,dz);

dCasA = decCasA;
rCasA = raCasA;
HCasA = (LST + h - rCasA).';
wCasA = @(dx,dy,dz) -cosd(dCasA)*sin(HCasA*pi/12)*X(dx) + sind(dCasA)*Y(dy,dz) + cosd(dCasA)*cos(HCasA*pi/12)*Z(dy,dz);
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
tGcygA = @(dx,dy,dz) reshape(wCygA(dx,dy,dz),1,tSize,6)/c;
matchedFilterCygA = @(dx,dy,dz,tI) sum(exp(-1i*2*pi*fc.*(tGcygA(dx,dy,dz)+tI)),3);
if antennaGainFactor
    % gainCygA = repmat(cosd(HCygA*15).',fSize,1);
    [~,gainCygA] = radiationPatternLofasm(fc,HCygA.');
    matchedFilterCygA = @(dx,dy,dz,tI) gainCygA .* matchedFilterCygA(dx,dy,dz,tI);
end

tGcasA = @(dx,dy,dz) reshape(wCasA(dx,dy,dz),1,tSize,6)/c;
matchedFilterCasA = @(dx,dy,dz,tI) sum(exp(-1i*2*pi*fc.*(tGcasA(dx,dy,dz)+tI)),3);
if antennaGainFactor
    % gainCasA = repmat(cosd(HCasA*15).',fSize,1);
    [~,gainCasA] = radiationPatternLofasm(fc,HCasA.');
    matchedFilterCasA = @(dx,dy,dz,tI) gainCasA .* matchedFilterCasA(dx,dy,dz,tI);
end

% Simulate AB Data
% ABSim = conj(matchedFilterCygA)+conj(matchedFilterCasA);
%---------------------------------------------------------------------%


%---------------------------------------------------------------------%
cygAPeakIntensity = @(dx,dy,dz,tI) sum(RlT.*matchedFilterCygA(dx,dy,dz,tI),2);
% figure; plot(fc,unwrap(angle(cygAPeakIntensity)));
casAPeakIntensity = @(dx,dy,dz,tI) sum(RlT.*matchedFilterCasA(dx,dy,dz,tI),2);
% figure; plot(fc,unwrap(angle(casAPeakIntensity)));


nFFT = 4*1024;
cygAPeakIntensityFFT = @(dx,dy,dz,tI) fftshift(fft(cygAPeakIntensity(dx,dy,dz,tI),nFFT));
casAPeakIntensityFFT = @(dx,dy,dz,tI) fftshift(fft(casAPeakIntensity(dx,dy,dz,tI),nFFT));


nIterations = 1;
dzEstMax = 10;
dxEstMax = 10;
dyEstMax = 10;
dzEst = 0;
dxEst = 0;
dyEst = 0;
tIEst = 0;

for nI = 1:nIterations
    %% Estimate dZ
    dzEst = linspace(dzEst-dzEstMax,dzEst+dzEstMax,50);
    cygAPeakIntensityFFTallXYZT = zeros(nFFT,length(dzEst));
    % casAPeakIntensityFFTallXYZT = zeros(nFFT,length(dzEst));
    for nZ = 1:length(dzEst)
        cygAPeakIntensityFFTallXYZT(:,nZ) = cygAPeakIntensityFFT(dxEst,dyEst,dzEst(nZ),tIEst);
        % casAPeakIntensityFFTallXYZT(:,nZ) = casAPeakIntensityFFT(dxEst,dyEst,dzEst(nZ),tIEst);
    end
    % mesh(dzEst,0:1023,abs(cygAPeakIntensityFFTallXYZT))
    % figure
    % mesh(dzEst,0:1023,abs(casAPeakIntensityFFTallXYZT))
    
    [~,dzInd] = max(max(abs(cygAPeakIntensityFFTallXYZT)));
    dzEst = dzEst(dzInd);
    dzEstMax = 3*dzEst;
    
    %% Estimate dX
    dxEst = linspace(dxEst-dxEstMax,dxEst+dxEstMax,50);
    cygAPeakIntensityFFTallXYZT = zeros(nFFT,length(dxEst));
    % casAPeakIntensityFFTallZ = zeros(nFFT,length(dzEst));
    for nX = 1:length(dxEst)
        cygAPeakIntensityFFTallXYZT(:,nX) = cygAPeakIntensityFFT(dxEst(nX),dyEst,dzEst,tIEst);
        % casAPeakIntensityFFTallZ(:,nZ) = casAPeakIntensityFFT(dxEst,dyEst,dzEst(nZ),tIEst);
    end
    [~,dxInd] = max(max(abs(cygAPeakIntensityFFTallXYZT)));
    dxEst = dxEst(dxInd);
    dxEstMax = 3*dxEst;
    
    %% Estimate dY
    dyEst = linspace(dyEst-dyEstMax,dyEst+dyEstMax,50);
    cygAPeakIntensityFFTallXYZT = zeros(nFFT,length(dyEst));
    % casAPeakIntensityFFTallZ = zeros(nFFT,length(dzEst));
    for nY = 1:length(dyEst)
        cygAPeakIntensityFFTallXYZT(:,nY) = cygAPeakIntensityFFT(dxEst,dyEst(nY),dzEst,tIEst);
        % casAPeakIntensityFFTallZ(:,nZ) = casAPeakIntensityFFT(dxEst,dyEst,dzEst(nZ),tIEst);
    end
    [~,dyInd] = max(max(abs(cygAPeakIntensityFFTallXYZT)));
    dyEst = dyEst(dyInd);
    dyEstMax = 3*dyEst;
end


%% Estimate Instrument Delay
[~,indTiCygA] = max(abs(cygAPeakIntensityFFTallXYZT));
% [~,indTiCasA] = max(abs(casAPeakIntensityFFTallZ));

tIEst = mean((indTiCygA-1-nFFT/2)/(nFFT*deltaF));
% tIEst = (indTiCasA-1-nFFT/2)/(nFFT*deltaF);

% figure; plot(indTiCygA)
% hold on
% plot(indTiCasA)


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


