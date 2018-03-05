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
% Outtriggger Antenna Coordinates
%-------------------------------------------------------------------------%
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

% Find the Local Siderial Time (LST) of Outtrigger Antenna
lstOut = 100.46 + 0.985647 * numOfDays + outLongMap + 15*startHour;
lstOut = mod(lstOut,360);
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
h = [raHl hT raHr] * pi/12; % in radians
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination
%-------------------------------------------------------------------------%
dec = linspace(0,80,50);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Delay Parameter
%-------------------------------------------------------------------------%
uDistance = 3; % Zenith baseline in meters
load('LoFASMAntenna60MHz.mat')
earthRad    = 6371008.7714; % equatorial radius (meters)
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gain Function of LoFASM Antenna
%---------------------------------------------------------------------%
gainFunction = zeros(length(dec),length(h));
%---------------------------------------------------------------------%

for nH = 1:length(h)
    
    [xOut,yOut,zOut] = sph2cart(mod(lstOut*pi/180+h(nH),2*pi),outLatMap*pi/180,earthRad+uDistance);
    
    rxLoFASMAntenna = phased.HeterogeneousConformalArray(...
        'ElementSet',{LoFASMAntenna},...
        'ElementIndices',1,...
        'ElementPosition',[xOut; yOut; zOut],...
        'ElementNormal',[wrapTo180(lstOut)+h(nH)*180/pi; outLatMap-90]); %  [azimuth;elevation]
%     figure; viewArray(rxLoFASMAntenna,'ShowIndex','All','ShowNormals',true)
%     pattern(rxLoFASMAntenna,60*1e6)
%     view(wrapTo180(raCygnusA*180/12)+90,decCygnusA)
%     drawnow();
    %---------------------------------------------------------------------%
    
    gainFunction(:,nH) = pattern(rxLoFASMAntenna,60*1e6,0,dec);
    
    release(rxLoFASMAntenna);
end
%-------------------------------------------------------------------------%