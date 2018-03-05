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
mainLatMap = 35.247227;
mainLongMap = -116.793272;

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

outLatMap = 35.247469;
outLongMap = -116.791509;
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
% Find the number of days from J2000.0 for 05:45:11 hrs UT on 2016 July 21th
%-------------------------------------------------------------------------%
startHour = 05 + 45/60 + 11/3600;
fractionDay = startHour / 24;
dayNum = 21;
daysSinceJanuary = 182;
daysSinceJ2000 = 5842.5;
numOfDays = daysSinceJ2000 + daysSinceJanuary + dayNum + fractionDay;

% Find the Local Siderial Time (LST) of Each Antenna
lstMain = 100.46 + 0.985647 * numOfDays + mainLongMap + 15*startHour;
lstMain = mod(lstMain,360);

lstMain1 = 100.46 + 0.985647 * numOfDays + mainLong1Map + 15*startHour;
lstMain1 = mod(lstMain1,360);
lstMain2 = 100.46 + 0.985647 * numOfDays + mainLong2Map + 15*startHour;
lstMain2 = mod(lstMain2,360);
lstMain3 = 100.46 + 0.985647 * numOfDays + mainLong3Map + 15*startHour;
lstMain3 = mod(lstMain3,360);
lstMain4 = 100.46 + 0.985647 * numOfDays + mainLong4Map + 15*startHour;
lstMain4 = mod(lstMain4,360);
lstMain5 = 100.46 + 0.985647 * numOfDays + mainLong5Map + 15*startHour;
lstMain5 = mod(lstMain5,360);
lstMain6 = 100.46 + 0.985647 * numOfDays + mainLong6Map + 15*startHour;
lstMain6 = mod(lstMain6,360);

lstOut = 100.46 + 0.985647 * numOfDays + outLongMap + 15*startHour;
lstOut = mod(lstOut,360);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Total Scanned Hour Duration
%-------------------------------------------------------------------------%
h = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
h = h * pi/12; % in radians
hSize = length(h);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination (in degrees)
%-------------------------------------------------------------------------%
dec = linspace(0,80,50);
dSize = length(dec);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Right Ascention (in degrees)
%-------------------------------------------------------------------------%
rA = linspace(-90,0,100);
rSize = length(rA);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Delay parameter
%-------------------------------------------------------------------------%
uDistance = 3; % Zenith baseline in meters
load('LoFASMAntenna60MHz.mat')
earthRad    = 6371008.7714; % equatorial radius (meters)
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Approximate Matched Filter with 1 Main Antenna is Used
%---------------------------------------------------------------------%
% tgApproximate = zeros(dSize,rSize,hSize);
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Exact Matched Filter with 6 Main Antenna is Used
%---------------------------------------------------------------------%
tgExact = zeros(dSize,rSize,6,hSize);
%---------------------------------------------------------------------%


for nH = 1:hSize
    for dD = 1:dSize       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % When Approximate Matched Filter with 1 Main Antenna is Used
        %---------------------------------------------------------------------%
        % [xMain,yMain,zMain] = sph2cart(mod(lstMain*pi/180+h(nH),2*pi),mainLatMap*pi/180,earthRad);
        %---------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % When Exact Matched Filter with 6 Main Antenna is Used
        %---------------------------------------------------------------------%
        [xMain1,yMain1,zMain1] = sph2cart(mod(lstMain1*pi/180+h(nH),2*pi),mainLat1Map*pi/180,earthRad);
        [xMain2,yMain2,zMain2] = sph2cart(mod(lstMain2*pi/180+h(nH),2*pi),mainLat2Map*pi/180,earthRad);
        [xMain3,yMain3,zMain3] = sph2cart(mod(lstMain3*pi/180+h(nH),2*pi),mainLat3Map*pi/180,earthRad);
        [xMain4,yMain4,zMain4] = sph2cart(mod(lstMain4*pi/180+h(nH),2*pi),mainLat4Map*pi/180,earthRad);
        [xMain5,yMain5,zMain5] = sph2cart(mod(lstMain5*pi/180+h(nH),2*pi),mainLat5Map*pi/180,earthRad);
        [xMain6,yMain6,zMain6] = sph2cart(mod(lstMain6*pi/180+h(nH),2*pi),mainLat6Map*pi/180,earthRad);
        %---------------------------------------------------------------------%
        
        [xOut,yOut,zOut] = sph2cart(mod(lstOut*pi/180+h(nH),2*pi),outLatMap*pi/180,earthRad+uDistance);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % When Approximate Matched Filter with 1 Main Antenna is Used
        %---------------------------------------------------------------------%
        %     rxLoFASM = phased.HeterogeneousConformalArray(...
        %         'ElementSet',{LoFASMAntenna},...
        %         'ElementIndices',[1 2],...
        %         'ElementPosition',[xMain xOut; yMain yOut; zMain zOut],...
        %         'ElementNormal',[wrapTo180([lstMain lstOut])+h(nH)*180/pi; [mainLatMap outLatMap]-90]); %  [azimuth;elevation]
        %     % viewArray(rxLoFASM,'ShowIndex','All','ShowNormals',true)
        %---------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % When Exact Matched Filter with 6 Main Antenna is Used
        %---------------------------------------------------------------------%
        rxLoFASM = phased.HeterogeneousConformalArray(...
            'ElementSet',{LoFASMAntenna},...
            'ElementIndices',[1 2 3 4 5 6 7],...
            'ElementPosition',[xMain1 xMain2 xMain3 xMain4 xMain5 xMain6 xOut;...
            yMain1 yMain2 yMain3 yMain4 yMain5 yMain6 yOut;...
            zMain1 zMain2 zMain3 zMain4 zMain5 zMain6 zOut],...
            'ElementNormal',[wrapTo180([lstMain1 lstMain2 lstMain3 lstMain4 lstMain5 lstMain6 lstOut])+h(nH)*180/pi;...
            [mainLat1Map mainLat2Map mainLat3Map mainLat4Map mainLat5Map mainLat6Map outLatMap]-90]); %  [azimuth;elevation]
%         viewArray(rxLoFASM,'ShowIndex','All','ShowNormals',true)
%         view(wrapTo180(raCygnusA*180/12)+90,decCygnusA)
%         drawnow();
        %---------------------------------------------------------------------%
        
        delayLoFASM = phased.ElementDelay('SensorArray',rxLoFASM);
        tau = delayLoFASM([rA;dec(dD)*ones(1,length(rA))]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % When Approximate Matched Filter with 1 Main Antenna is Used
        %---------------------------------------------------------------------%
        %     tgApproximate(dD,:,nH) = diff(tau);
        %---------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % When Exact Matched Filter with 6 Main Antenna is Used
        %---------------------------------------------------------------------%
        for nA = 1:6
            tgExact(dD,:,nA,nH) = tau(7,:)-tau(nA,:);
        end
        %---------------------------------------------------------------------%
        
        release(delayLoFASM);
        
    end
end