tStart = 1;
fStart = 550; % 140; % 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz
fStop = 550; %800; % 1024 = 100 MHz, 820 = 80 MHz, 417 = 40 MHz

SimulationMode = false;

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
    ABlTerm = ABSim;
    %-------------------------------------------------------------------------%
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Real Data
    %-------------------------------------------------------------------------%
    if nSAv == 447
        load('ABlTermPZ');
    else
        load('AB_Padded.mat')
        ABlTerm = longTermAverage(AB,nSAv);      
    end
end

ABlTerm = whiteningApproach1(ABlTerm);
ABlT = ABlTerm;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract Related Frequency and Time Range
%-------------------------------------------------------------------------%
ABlT = ABlT(fStart:fStop,tStart:end);
[fSize, tSize] = size(ABlT);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna Coordinates
%-------------------------------------------------------------------------%
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

% lstMain1 = mod(mainLong1Map + GST*15,360);
% lstMain2 = mod(mainLong2Map + GST*15,360);
% lstMain3 = mod(mainLong3Map + GST*15,360);
% lstMain4 = mod(mainLong4Map + GST*15,360);
% lstMain5 = mod(mainLong5Map + GST*15,360);
% lstMain6 = mod(mainLong6Map + GST*15,360);

LST = mod(outLongMap/15 + GST,24);
h = (0 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
% LST = LST + h;
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

% innerNED = [2.205*[-1, -2, -1, 1, 2, 1];...
%             2.205*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
%             zeros(1,6)];

innerNED = [2.205*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
            2.205*[-1, -2, -1, 1, 2, 1];...
            zeros(1,6)];

% plot(innerNED(2,:),innerNED(1,:),'.')
        
% innerNED2 = innerNED;
% innerNED2(2,:) = innerNED(1,:);
% innerNED2(1,:) = innerNED(2,:);
% innerNED = innerNED2;

innerNED = rotNed2Body * innerNED + [northMO;eastMO;downMO];

x = innerNED(2,:);
y = innerNED(1,:);
z = -innerNED(3,:);

% plot(x,y,'o') x: towards the east, y: toward the north
% sum((innerNED(:,1)-innerNED(:,2)).^2)^0.5
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Scanned Right Ascention
%-------------------------------------------------------------------------%
rA = linspace(18,24,200); % in hours
rSize = length(rA);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination
%-------------------------------------------------------------------------%
% dec = linspace(0,80,200).'; % in degrees
dec = decCygnusA;
dSize = length(dec);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
fc = (0.09765625*(fStart-1):0.09765625:0.09765625*(fStop-1))'*1e6;
c = physconst('LightSpeed');
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baselines
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = x; X = reshape(X,1,1,6);
Y = z*sind(l) + y*cosd(l); Y = reshape(Y,1,1,6);
Z = z*cosd(l) - y*sind(l); Z = reshape(Z,1,1,6);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare AB for IFFT
%-------------------------------------------------------------------------%
ABlT = reshape(transpose(ABlT),1,tSize,fSize);
fc = reshape(fc,1,1,fSize);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Instrument Delay
%-------------------------------------------------------------------------%
if SimulationMode
    tI = 0;
else
    nsF = 67; % 16; 68 ; 135; 69.5; % Best Value is 67
    tI = nsF/(100e6)*ones(1,6);
%     tI = ones(1,6)*6.4000e-07;
%     tI = 0.0068; % T_newton; % ones(1,6)*6.7000e-07; % 4.0371e-07; %  % 0.0068; % -0.0011; % 
end
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare Intensity
%-------------------------------------------------------------------------%
intensityMF = zeros(dSize,rSize);
%-------------------------------------------------------------------------%


if (dSize == 1)
    figure
    xlabel('RA (Hour)')
end

for nH = 1:length(h)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate Theoric Delay
    %-------------------------------------------------------------------------%
    % w = zeros(length(dec),length(H),6);
    H = LST + h(nH) - rA;
    w = -cosd(dec)*sin(H*pi/12).*X + sind(dec).*Y + cosd(dec)*cos(H*pi/12).*Z;
    %-------------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% When Exact Matched Filter with 6 Main Antenna is Used
    %---------------------------------------------------------------------%
    matchedFilter = 0;
    for nA = 1:6
        tG = squeeze(w(:,:,nA))/c;
        matchedFilter = matchedFilter + exp(-2*pi*1i*fc.*(tG-tI(nA)));
    end
    %---------------------------------------------------------------------%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 2D FFT of AB Data
    ABlTfft = fft2(ABlT(:,nH,:),dSize,rSize);
    %-------------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate Intensity
    %-------------------------------------------------------------------------%
    intensityMF = intensityMF + sum(ifft2(ABlTfft.*fft2(matchedFilter)),3);
    %-------------------------------------------------------------------------%
    if (dSize == 1)
        plot(rA,abs(intensityMF))
        drawnow;
    end
    
end

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
