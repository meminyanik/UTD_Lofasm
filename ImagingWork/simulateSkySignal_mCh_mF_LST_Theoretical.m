addpath(genpath('C:\D\UTD\Research\LoFASM\harddiskRecordedData'));
addpath(genpath('C:\D\UTD\Research\LoFASM\ImagingWork'))

%-------------------------------------------------------------------------%
% Initial Variables
%-------------------------------------------------------------------------%
fRange = 1:1024;
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

antennaGainFactor = false;
SourceModelMode = 3; % 1: CygA, 2: Cas A, 3: CygA+CasA

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
% Load AB Long term 
%-------------------------------------------------------------------------%
if nSAv == 447
    load('ABlTermPZ');
else
    load('AB_Padded.mat')
    ABlTerm = longTermAverage(AB,nSAv);
end
ABlTerm = whiteningApproach1(ABlTerm);
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
%%% Calculate Baselines and Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = x;
Y = z*sind(l) + y*cosd(l);
Z = z*cosd(l) - y*sind(l);

dCygA = decCygnusA;
rCygA = raCygnusA;
HCygA = (LST + h - rCygA).';
wCygA = -cosd(dCygA)*sin(HCygA*pi/12)*X + sind(dCygA)*Y + cosd(dCygA)*cos(HCygA*pi/12)*Z;

dCasA = decCasA;
rCasA = raCasA;
HCasA = (LST + h - rCasA).';
wCasA = -cosd(dCasA)*sin(HCasA*pi/12)*X + sind(dCasA)*Y + cosd(dCasA)*cos(HCasA*pi/12)*Z;
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
tI = repmat(-6.9000e-07,1,6); % -6.9000e-07;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate AB Sim
%---------------------------------------------------------------------%
ABSim = zeros(fSize,tSize);
ABSimEachAntenna = zeros(fSize,tSize);
if antennaGainFactor
        gainCygA = 8100*cosd(HCygA*15).'; % ones(1,length(HCygA));
        gainCasA = 11000*cosd(HCasA*15).'; % ones(1,length(HCasA));
else
    gainCygA = 1;
    gainCasA = 1;
end

antennaGain = [2, 3, 5, 1, 1, 1];
% antennaGain = [1, 1, 1, 1, 1, 1];
antennaGain = antennaGain/norm(antennaGain);

%%
for nA = 1:6
    tGCygA = (wCygA(:,nA)/c).';
    tGCasA = (wCasA(:,nA)/c).';
    
    switch SourceModelMode
        case 1
            ABSimEachAntenna = gainCygA.^2.*exp(1i*2*pi*fc.'*(tGCygA+tI(nA)));
        case 2
            ABSimEachAntenna = gainCasA.^2.*exp(1i*2*pi*fc.'*(tGCasA+tI(nA)));
        case 3
            ABSimEachAntenna = gainCygA.^2.*exp(1i*2*pi*fc.'*(tGCygA+tI(nA))) + gainCasA.^2.*exp(1i*2*pi*fc.'*(tGCasA+tI(nA)));
    end
    
    ABSim = ABSim + antennaGain(nA)*ABSimEachAntenna;
end
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Further Analysis
%-------------------------------------------------------------------------%
ABlTifft = ifft(ABlT);
ABSimifft = ifft(ABSim);
mesh(abs(ABlTifft(20:200,:)))
figure
mesh(abs(ABSimifft(20:200,:)))

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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Old Version
% %-------------------------------------------------------------------------%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Set Antenna Parameters
% %-------------------------------------------------------------------------%
% load('LoFASMAntenna60MHz.mat')
% %-------------------------------------------------------------------------%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Set Transmitted Signal
% %-------------------------------------------------------------------------%
% % Create transmit antenna
% txAntenna = phased.IsotropicAntennaElement('FrequencyRange',[10,80]*1e6);
% 
% % Create transmitting signal
% % f1 = 600.0;
% nS = 2048*8192; % Number of samples in each timing period
% % fS = 200e6; % Sampling Rate
% % t = linspace(0,0.083886,nS);
% % xS = cos(2*pi*t*f1)';
% xS = randn(1,nS)' + 1i*randn(1,nS)';
% %-------------------------------------------------------------------------%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Calculate Simulated Signal
% %-------------------------------------------------------------------------%
% for hH = 1:length(h)
%     for fF = 1:length(fc)
%         
%         % Create radiator
%         radiator = phased.Radiator('Sensor',txAntenna,'OperatingFrequency',fc(fF));
%         xT = radiator(xS,[-180;0]);
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % When Exact Matched Filter with 6 Main Antenna is Used
%         %---------------------------------------------------------------------%
%         [xMain1,yMain1,zMain1] = sph2cart(mod(lstMain1*pi/180+h(hH),2*pi),mainLat1Map*pi/180,earthRad);
%         [xMain2,yMain2,zMain2] = sph2cart(mod(lstMain2*pi/180+h(hH),2*pi),mainLat2Map*pi/180,earthRad);
%         [xMain3,yMain3,zMain3] = sph2cart(mod(lstMain3*pi/180+h(hH),2*pi),mainLat3Map*pi/180,earthRad);
%         [xMain4,yMain4,zMain4] = sph2cart(mod(lstMain4*pi/180+h(hH),2*pi),mainLat4Map*pi/180,earthRad);
%         [xMain5,yMain5,zMain5] = sph2cart(mod(lstMain5*pi/180+h(hH),2*pi),mainLat5Map*pi/180,earthRad);
%         [xMain6,yMain6,zMain6] = sph2cart(mod(lstMain6*pi/180+h(hH),2*pi),mainLat6Map*pi/180,earthRad);
%         %---------------------------------------------------------------------%
%         
%         [xOut,yOut,zOut] = sph2cart(mod(lstOut*pi/180+h(hH),2*pi),outLatMap*pi/180,earthRad+uDistance);
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % When Exact Matched Filter with 6 Main Antenna is Used
%         %---------------------------------------------------------------------%
%         rxLoFASM = phased.HeterogeneousConformalArray(...
%             'ElementSet',{LoFASMAntenna},...
%             'ElementIndices',[1 1 1 1 1 1 1],...
%             'ElementPosition',[xMain1 xMain2 xMain3 xMain4 xMain5 xMain6 xOut;...
%             yMain1 yMain2 yMain3 yMain4 yMain5 yMain6 yOut;...
%             zMain1 zMain2 zMain3 zMain4 zMain5 zMain6 zOut],...
%             'ElementNormal',[wrapTo180([lstMain1 lstMain2 lstMain3 lstMain4 lstMain5 lstMain6 lstOut])+h(hH)*180/pi;...
%             [mainLat1Map mainLat2Map mainLat3Map mainLat4Map mainLat5Map mainLat6Map outLatMap]-90]); %  [azimuth;elevation]
%         % viewArray(rxLoFASM,'ShowIndex','All','ShowNormals',true)
%         %---------------------------------------------------------------------%
%         
%                 
%         collector = phased.Collector('Sensor',rxLoFASM,'OperatingFrequency',fc(fF));
%         yR = collector(xT,anglesSimSource);
%         
%         aaR = sum(yR(:,1:6),2);
%         bbR = yR(:,7);
%         aaR = reshape(aaR,2048,8192);
%         bbR = reshape(bbR,2048,8192);
%         
%         aaRfft = fft(aaR);
%         bbRfft = fft(bbR);
%         
%         aaR = aaRfft .* conj(aaRfft);
%         aaR = sum(aaR,2)/8192;
%         aaR = aaR(1:1024);
%         
%         abR = aaRfft .* conj(bbRfft);
%         abR = sum(abR,2)/8192;
%         abR = abR(1:1024);
%         
%         bbR = bbRfft .* conj(bbRfft);
%         bbR = sum(bbR,2)/8192;
%         bbR = bbR(1:1024);
%         
%         AASim(:,hH) = aaR;
%         BBSim(:,hH) = bbR;
%         ABSim(:,hH) = abR;
%         
%         release(radiator);
%         release(collector);
%     end
% end