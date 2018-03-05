%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the number of days from J2000.0 for 05:45:11 hrs UT on 2016 July 21th
%-------------------------------------------------------------------------%
startHour = 05 + 45/60 + 11/3600;
fractionDay = startHour / 24;
dayNum = 21;
daysSinceJanuary = 182;
daysSinceJ2000 = 5842.5;
numOfDays = daysSinceJ2000 + daysSinceJanuary + dayNum + fractionDay;
%-------------------------------------------------------------------------%

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
% RA of Sky Sources
%-------------------------------------------------------------------------%
raCygnusA = 19 + 59/60 + 28.3/3600; % in Hours
decCygnusA = 40 + 44/60 + 02/3600; % in Degrees
% Cygnus A: 19h59m28.3s +40d44m02s
raCasA = 23 + 23/60 + 27.9/3600; % in Hours
decCasA = 58 + 48/60 + 42/3600; % in Degrees
% Cas A: 23h23m27.9s +58d48m42s
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LST of the station (outtrigger) at the beginning of record
%-------------------------------------------------------------------------%
LST = 100.46 + 0.985647 * numOfDays + outLongMap + 15*startHour;
LST = mod(LST,360);

% LST of the station during the whole record
tSize = 504;
hRecording = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
LSTRecording = LST + 15*hRecording;
LSTRecording = mod(LSTRecording,360);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HA of sky sources during the whole record
%-------------------------------------------------------------------------%
haCygnusA = LSTRecording - raCygnusA;
haCygnusA = mod(haCygnusA,360);
haCasA = LSTRecording - raCasA;
haCasA = mod(haCasA,360);
%-------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Antenna Coordinates
%-------------------------------------------------------------------------%
wgs84 = wgs84Ellipsoid('meters');
[main1North,main1East,main1Down] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[main2North,main2East,main2Down] = geodetic2ned(mainLat2Map,mainLong2Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[main3North,main3East,main3Down] = geodetic2ned(mainLat3Map,mainLong3Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[main4North,main4East,main4Down] = geodetic2ned(mainLat4Map,mainLong4Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[main5North,main5East,main5Down] = geodetic2ned(mainLat5Map,mainLong5Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
[main6North,main6East,main6Down] = geodetic2ned(mainLat6Map,mainLong6Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);

% Baseline: x-towards the east, y-toward the north, and z-toward the zenith
localXYZ = [main1East, main1North, -main1Down;...
            main2East, main2North, -main2Down;...
            main3East, main3North, -main3Down;...
            main4East, main4North, -main4Down;...
            main5East, main5North, -main5Down;...
            main6East, main6North, -main6Down];
localXYZ = localXYZ';

% Intermediate transformation
intermediateTransform = [1, 0, 0;...
                0, cosd(outLatMap), sind(outLatMap);...
                0, -sind(outLatMap), cosd(outLatMap)];
intermediateXYZ = intermediateTransform * localXYZ;
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w-Plane for CygnusA
%-------------------------------------------------------------------------%
wPlaneTransform = [-cosd(decCygnusA)*sind(haCygnusA); ...
                    sind(decCygnusA)*ones(1,length(haCygnusA));...
                    cosd(decCygnusA)*cosd(haCygnusA)];

wPlane = wPlaneTransform' * intermediateXYZ;

% Time delay
c = physconst('LightSpeed');
timeDelay = wPlane' / c;
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated Sky Signal
%-------------------------------------------------------------------------%
fStart = 1; % 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz
fStop = 1024; % 1024 = 100 MHz, 820 = 80 MHz, 417 = 40 MHz
fc = (0.09765625*(fStart-1):0.09765625:0.09765625*(fStop-1))'*1e6;

skySignal = 0;
for nA = 1:6
    skySignal = skySignal + exp(2*pi*1i*fc.*timeDelay(nA,:));
end
%-------------------------------------------------------------------------%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
% nSAv = 50;
%-------------------------------------------------------------------------%

recordStart = (5 + 45/60 + 11/3600);
% recordStart = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For Long Term Correlation
%-------------------------------------------------------------------------%
hT = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
raLd = -20; raRd = 50;
raLs = round(raLd*12/180*3600/(nSAv*0.083886));
raRs = round(raRd*12/180*3600/(nSAv*0.083886));

raHl = (raLs*nSAv*0.083886 : nSAv*0.083886 : -1*nSAv*0.083886)/3600;
raHr = (tSize*nSAv*0.083886 : nSAv*0.083886 : (tSize+raRs)*nSAv*0.083886)/3600;
h = [raHl hT raHr];
h = (h + recordStart) * pi/12;
% h = flip(h);
%-------------------------------------------------------------------------%

dec = linspace(0,80,50)';

c = physconst('LightSpeed');

tgExactM = zeros(length(dec),6,length(h));
for nA = 1:6
    [arclen,az] = distance(mainLatMap(nA),mainLongMap(nA),outLatMap,outLongMap,referenceEllipsoid('earth','m'));
    eDistance = arclen*sind(az); % East-West baseline in meters
    nDistance = arclen*cosd(az); % North-South baseline in meters
    uDistance = 3; % Zenith baseline in meters
    
    
    enu = [eDistance;nDistance;uDistance];
    rotM = [0 -sin(outLatMap) cos(outLatMap);1 0 0;0 cos(outLatMap) sin(outLatMap)];
    Lxyz = rotM * enu;
    
    tgExactM(:,nA,:) = 1/c * (Lxyz(1)*cosd(dec)*cos(h) - Lxyz(2)*cosd(dec)*sin(h) + Lxyz(3)*sind(dec)*ones(1,length(h)));
end