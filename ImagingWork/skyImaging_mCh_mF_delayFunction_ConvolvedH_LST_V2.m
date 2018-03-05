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
h = [raHl hT raHr] * 180/12; % in degrees
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination
%-------------------------------------------------------------------------%
dec = linspace(0,80,50);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Parameters
%-------------------------------------------------------------------------%
load('LoFASMAntenna60MHz.mat')
wgs84 = wgs84Ellipsoid('meters');
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Exact Matched Filter with 6 Main Antenna is Used
%---------------------------------------------------------------------%
tgExact = zeros(length(dec),6,length(h));
%---------------------------------------------------------------------%

[xMain1,yMain1,zMain1] = geodetic2ecef(wgs84,mainLat1Map,mainLong1Map,mainAltMap);
[xMain2,yMain2,zMain2] = geodetic2ecef(wgs84,mainLat2Map,mainLong2Map,mainAltMap);
[xMain3,yMain3,zMain3] = geodetic2ecef(wgs84,mainLat3Map,mainLong3Map,mainAltMap);
[xMain4,yMain4,zMain4] = geodetic2ecef(wgs84,mainLat4Map,mainLong4Map,mainAltMap);
[xMain5,yMain5,zMain5] = geodetic2ecef(wgs84,mainLat5Map,mainLong5Map,mainAltMap);
[xMain6,yMain6,zMain6] = geodetic2ecef(wgs84,mainLat6Map,mainLong6Map,mainAltMap);

[xOut,yOut,zOut] = geodetic2ecef(wgs84,outLatMap,outLongMap,outAltMap);

xyzMain1Initial = rotz(GST*15)*[xMain1;yMain1;zMain1];
xyzMain2Initial = rotz(GST*15)*[xMain2;yMain2;zMain2];
xyzMain3Initial = rotz(GST*15)*[xMain3;yMain3;zMain3];
xyzMain4Initial = rotz(GST*15)*[xMain4;yMain4;zMain4];
xyzMain5Initial = rotz(GST*15)*[xMain5;yMain5;zMain5];
xyzMain6Initial = rotz(GST*15)*[xMain6;yMain6;zMain6];

xyzOutInitial = rotz(GST*15)*[xOut;yOut;zOut];

for nH = 1:length(h)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % When Exact Matched Filter with 6 Main Antenna is Used
    %---------------------------------------------------------------------%
    xyzMain1 = rotz(h(nH)) * xyzMain1Initial;
    xyzMain2 = rotz(h(nH)) * xyzMain2Initial;
    xyzMain3 = rotz(h(nH)) * xyzMain3Initial;
    xyzMain4 = rotz(h(nH)) * xyzMain4Initial;
    xyzMain5 = rotz(h(nH)) * xyzMain5Initial;
    xyzMain6 = rotz(h(nH)) * xyzMain6Initial;
    
    xyzOut = rotz(h(nH)) * xyzOutInitial;
    %---------------------------------------------------------------------%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % When Exact Matched Filter with 6 Main Antenna is Used
    %---------------------------------------------------------------------%
    rxLoFASM = phased.HeterogeneousConformalArray(...
        'ElementSet',{LoFASMAntenna},...
        'ElementIndices',[1 2 3 4 5 6 7],...
        'ElementPosition',[xyzMain1(1) xyzMain2(1) xyzMain3(1) xyzMain4(1) xyzMain5(1) xyzMain6(1) xyzOut(1);...
        xyzMain1(2) xyzMain2(2) xyzMain3(2) xyzMain4(2) xyzMain5(2) xyzMain6(2) xyzOut(2);...
        xyzMain1(3) xyzMain2(3) xyzMain3(3) xyzMain4(3) xyzMain5(3) xyzMain6(3) xyzOut(3)],...
        'ElementNormal',[wrapTo180([lstMain1 lstMain2 lstMain3 lstMain4 lstMain5 lstMain6 lstOut]+h(nH));...
        [mainLat1Map mainLat2Map mainLat3Map mainLat4Map mainLat5Map mainLat6Map outLatMap]-90]); %  [azimuth;elevation]
%     viewArray(rxLoFASM,'ShowIndex','All','ShowNormals',true)
%     view(wrapTo180(raCygnusA*180/12)+90,decCygnusA)
%     drawnow();
    %---------------------------------------------------------------------%
    
    delayLoFASM = phased.ElementDelay('SensorArray',rxLoFASM);
    tau = delayLoFASM([wrapTo180(lstOut+raRd)*ones(1,length(dec));dec]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         % When Approximate Matched Filter with 1 Main Antenna is Used
    %---------------------------------------------------------------------%
    %     tgApproximate(:,nH) = diff(tau);
    %---------------------------------------------------------------------%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         % When Exact Matched Filter with 6 Main Antenna is Used
    %---------------------------------------------------------------------%
    for nA = 1:6
        tgExact(:,nA,nH) = tau(7,:)-tau(nA,:);
    end
    %---------------------------------------------------------------------%
    
    release(delayLoFASM);
end
%-------------------------------------------------------------------------%