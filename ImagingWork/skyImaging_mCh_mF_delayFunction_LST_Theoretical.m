%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
% nSAv = 50;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load AB Long Term 
%-------------------------------------------------------------------------%
if nSAv == 447
    load('ABlTermPZ.mat')
else
    load('AB_Padded.mat')
    ABlTerm = longTermAverage(AB,nSAv);
end

ABlT = ABlTerm;
[~, tSize] = size(ABlT);
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
LST = LST + h;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baseline from Map Locations
%-------------------------------------------------------------------------%
wgs84 = wgs84Ellipsoid('meters');

arrayCenterLat = (mainLat2Map+mainLat5Map)/2;
arrayCenterLon = (mainLong2Map+mainLong5Map)/2;

[northMO,eastMO,downMO] = geodetic2ned(arrayCenterLat,arrayCenterLon,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);

innerRingRotation = azimuth(mainLat2Map,mainLong2Map,mainLat5Map,mainLong5Map,wgs84); % Counter-clockwise - Degrees
rotNed2Body = angle2dcm( innerRingRotation*pi/180, 0, 0 );

innerNED = [2.205*[-1, -2, -1, 1, 2, 1];...
            2.205*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
            zeros(1,6)];

innerNED = rotNed2Body * innerNED + [northMO;eastMO;downMO];

x = innerNED(2,:);
y = innerNED(1,:);
z = innerNED(3,:);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Right Ascention
%-------------------------------------------------------------------------%
RA = linspace(18,24,300); % in hours
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination
%-------------------------------------------------------------------------%
DEC = linspace(30,70,200);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = x;
Y = z*sind(l) + y*cosd(l);
Z = z*cosd(l) - y*sind(l);

w = zeros(length(DEC),length(RA),6,length(h));
for dD = 1:length(DEC)
    d = DEC(dD);
    for rR = 1:length(RA)
        H = ((LST - RA(rR))*pi/12).';
        w(dD,rR,:,:) = (-cosd(d)*sin(H)*X + sind(d)*Y + cosd(d)*cos(H)*Z).';
    end
end
%-------------------------------------------------------------------------%

c = physconst('LightSpeed');
tgExact = w/c;