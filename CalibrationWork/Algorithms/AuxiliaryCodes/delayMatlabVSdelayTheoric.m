%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load Matlab Delay
%-------------------------------------------------------------------------%
load('tgExact_LST_V2.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination (in degrees)
%-------------------------------------------------------------------------%
dec = linspace(0,80,50);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Right Ascention (in degrees)
%-------------------------------------------------------------------------%
rA = linspace(-90,0,100);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Exact Matched Filter with 6 Main Antenna is Used
%---------------------------------------------------------------------%
[dSize,rSize,nAntenna,hSize] = size(tgExact);
%---------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna Coordinates
%-------------------------------------------------------------------------%
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
mainAltMap = distdim(3547,'ft','m');

outLatMap = 35.247469;
outLongMap = -116.791509;
uDistance = 3; % Zenith baseline in meters
outAltMap = distdim(3525,'ft','m') + uDistance;
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
lstTurn = 100.46 + 0.985647 * numOfDays + 15*startHour;
lstTurn = mod(lstTurn,360);

% Find the Local Siderial Time (LST) of Each Antenna
lstMain1 = mod(mainLong1Map + lstTurn,360);
lstMain2 = mod(mainLong2Map + lstTurn,360);
lstMain3 = mod(mainLong3Map + lstTurn,360);
lstMain4 = mod(mainLong4Map + lstTurn,360);
lstMain5 = mod(mainLong5Map + lstTurn,360);
lstMain6 = mod(mainLong6Map + lstTurn,360);

lstOut = mod(outLongMap + lstTurn,360);
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
% Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Total Scanned Hour Duration
%-------------------------------------------------------------------------%
tSize = 504;
h = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
h = h * 180/12; % in degrees
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Theoric Delay
%-------------------------------------------------------------------------%
wgs84 = wgs84Ellipsoid('meters');
% Outtrigger is the reference
[north1,east1,down1] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[north2,east2,down2] = geodetic2ned(mainLat2Map,mainLong2Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[north3,east3,down3] = geodetic2ned(mainLat3Map,mainLong3Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[north4,east4,down4] = geodetic2ned(mainLat4Map,mainLong4Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[north5,east5,down5] = geodetic2ned(mainLat5Map,mainLong5Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[north6,east6,down6] = geodetic2ned(mainLat6Map,mainLong6Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);

l = outLatMap;
L = outLongMap;

x = [east1, east2, east3, east4, east5, east6];
y = [north1, north2, north3, north4, north5, north6];
z = -1*[down1, down2, down3, down4, down5, down6];

X = x;
Y = z*sind(l) + y*cosd(l);
Z = z*cosd(l) - y*sind(l);

dIndx = 26;
rIndx = 34;
d = decCygnusA; % dec(dIndx);
r = raCygnusA * 15; % rA(rIndx);
H = (lstOut + h - r).';

w = -cosd(d)*sind(H)*X + sind(d)*Y + cosd(d)*cosd(H)*Z;

% c = physconst('LightSpeed');
% plot(w(:,2)/c)
% hold on
% plot(squeeze(tgExact(dIndx,rIndx,2,:)))


%% Simulate AB Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
f = (0.09765625*0:0.09765625:0.09765625*1023)'*1e6;
c = physconst('LightSpeed');
%-------------------------------------------------------------------------%

w = reshape(w,1,tSize,6);
ABSim = sum(exp(1i*2*pi*f.*w/c),3);

plot(real(ABSim(550,:)))