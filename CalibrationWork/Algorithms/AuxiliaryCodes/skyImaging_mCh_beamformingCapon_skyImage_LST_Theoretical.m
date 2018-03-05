addpath(genpath('C:\D\UTD\Research\LoFASM\CalibrationWork'))
addpath(genpath('C:\D\UTD\Research\LoFASM\harddiskRecordedData'));
addpath(genpath('C:\D\UTD\Research\LoFASM\ImagingWork'));

%-------------------------------------------------------------------------%
% Initial Variables
%-------------------------------------------------------------------------%
% fRange = [300:375 450:500 550:650];
fRange = [300:800];
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

SimulationMode = false;
Whitening = true;

useCalibrationPerFrequency = false;

StationID = 4;

%-------------------------------------------------------------------------%
% Data Record Time
%-------------------------------------------------------------------------%
% Data Record Start Time [h,m,s,d,m,y]
dataRecordTime = [5,45,11,21,7,2016];

% dataRecordTime = [20,35,23,4,12,2017];
% dataRecordTime = [19,51,50,11,11,2017];
% dataRecordTime = [23,53,49,11,11,2017];

% dataRecordTime = [16,59,59,26,4,2017]; % start time
% dataRecordTime = [11,28,41,27,4,2017]; % ~18:00 LST
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
% nSAv = 50;
% nSAv = 1;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Real Data
%-------------------------------------------------------------------------%
if (~SimulationMode)
    if nSAv == 447
        load('AAlTermPZ');
        load('ABlTermPZ');
        load('BBlTermPZ');
    else
        load('AA_Padded.mat')
        load('AB_Padded.mat')
        load('BB_Padded.mat')
        AAlTerm = longTermAverage(AA,nSAv);
        ABlTerm = longTermAverage(AB,nSAv);
        BBlTerm = longTermAverage(BB,nSAv);
    end
    
    if Whitening
        ABlTerm = whiteningApproach1(ABlTerm);
    end
    
%     ABlTerm((abs(ABlTerm)./(AAlTerm+BBlTerm)) > 0.1) = 0;
%     AAlTerm((abs(ABlTerm)./(AAlTerm+BBlTerm)) > 0.1) = 0;
%     BBlTerm((abs(ABlTerm)./(AAlTerm+BBlTerm)) > 0.1) = 0;
    
    AAlT = AAlTerm; % clear AAlTerm
    ABlT = ABlTerm; % clear ABlTerm
    BBlT = BBlTerm; % clear BBlTerm

else
    AAlT = AASim;
    ABlT = ABSim;
    BBlT = BBSim;
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract Related Frequency and Time Range
%-------------------------------------------------------------------------%
AAlT = AAlT(fRange,:);
ABlT = ABlT(fRange,:);
BBlT = BBlT(fRange,:);
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
    
elseif StationID == 2
    mainLat2Map = 34.0787;
    mainLong2Map = -107.6183;
    
    mainLat5Map = 0;
    mainLong5Map = 0;
    mainAltMap = 0;
    
    outLatMap = 0;
    outLongMap = 0;
    outAltMap = 0;
    
elseif StationID == 3   %%%  Green Bank
    mainLat2Map = 38.429034;
    mainLong2Map = -79.846285;
    
    mainLat5Map = 38.429060;
    mainLong5Map = -79.846189;
    mainAltMap = 807;
    
    outLatMap = 38.429791;
    outLongMap = -79.843287;
    outAltMap = 804;
    
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
raCygnusA = 19 + 59/60 + 28.3/3600; % in Hours (19.9912)
decCygnusA = 40 + 44/60 + 02/3600; % in Degrees (40.7339)
% Cygnus A: 19h59m28.3s +40d44m02s
raCasA = 23 + 23/60 + 27.9/3600; % in Hours (23.3911)
decCasA = 58 + 48/60 + 42/3600; % in Degrees (58.8117)
% Cas A: 23h23m27.9s +58d48m42s
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the LST
%-------------------------------------------------------------------------%
% LST = calculate_LST(longitude,hour,minute,second,day,month,year);
dataRecordTime = num2cell(dataRecordTime);
LSTest = calculate_LST(outLongMap,dataRecordTime{:});
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

% plot(innerNED(2,:),innerNED(1,:),'o')

innerNED = rotNed2Body * innerNED + [northMO;eastMO;downMO];

dxEst = 0.3030; % 0.2041; for (LoFASM1) % or -0.2041; % 0.3030 for (LoFASM4)
dyEst = 0.6566; % 0.6122; for (LoFASM1) % or 0.2041; % 0.6566 for (LoFASM4)
dzEst = 2.1212; % -0.2041; for (LoFASM1) % or 1.8367; % 2.1212 for (LoFASM4)

x = innerNED(2,:) + dxEst; % LoFASMIV: 
y = innerNED(1,:) + dyEst;
z = -innerNED(3,:) + dzEst; % LoFASMIV: 

% D = sqrt(x.^2+y.^2+z.^2);
% az = 270-atand(y./x);
% el = atand(z./D);

% [X,Y,Z] = @(l,L) ned2ecef(xNorth,yEast,zDown,l,L,h0,spheroid);

% plot(x,y,'o') x: towards the east, y: toward the north
% sum((innerNED(:,1)-innerNED(:,2)).^2)^0.5
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Scanned Right Ascention and Total Scanned Hour Angle
%-------------------------------------------------------------------------%
lSearch = linspace(20,50,50); % in degrees
rA = raCygnusA;
% rA = raCasA;
lSize = length(lSearch);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination
%-------------------------------------------------------------------------%
LSearch = linspace(15,20,50); % in hours
dec = decCygnusA;
% dec = decCasA;
LSize = length(LSearch);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baselines and Theoric Delay
%-------------------------------------------------------------------------%
% l = outLatMap;
% L = outLongMap;

X = x; X = reshape(X,1,1,6);
Y = @(l) z*sind(l) + y*cosd(l);
Z = @(l) z*cosd(l) - y*sind(l);

dY = @(l) z*cosd(l) - y*sind(l);
dZ = @(l) -z*sind(l) - y*cosd(l);




        

        GST = calculate_GST(dataRecordTime{:});
        LST = @(L) mod(L/15+GST,24);

        dCygA = decCygnusA;
        rCygA = raCygnusA;
        % HCygA = @(L) LST(L) + h - rCygA;
        % dHCygA = 1/15;
        HCygA = @(L) (LST(L) + h - rCygA)*pi/12; % L is in degrees
        dHCygA = 1;
        
        wCygA = @(l,L) -cosd(dCygA)*sin(HCygA(L))*X + sind(dCygA)*Y(l) + cosd(dCygA)*cos(HCygA(L))*Z(l);
        dwCygA_l = @(l,L) sind(dCygA)*dY(l) + cosd(dCygA)*cos(HCygA(L))*dZ(l);
        dwCygA_L = @(l,L) -cosd(dCygA)*cos(HCygA(L))*X - cosd(dCygA)*sin(HCygA(L))*Z(l);


        

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Instrument Delay
%-------------------------------------------------------------------------%
% tI = 67/(100e6);
tI = -6.9000e-07; % tIEst;
tI = repmat(tI,1,1,6);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize Sky Intensity
%-------------------------------------------------------------------------%
% intensityBF = zeros(dSize,rSize);
intensityMF = zeros(lSize,LSize);
% intensityMF_EachFreqTime = zeros(dSize,rSize,fSize,tSize);
% intensityCAPON = zeros(dSize,rSize);
% intensityCAPONAppr1 = zeros(dSize,rSize); % Capon with summation accross frequency after inversion
% intensityCAPONAppr2 = zeros(dSize,rSize); % Capon with summation accross frequency before inversion
% intensityCAPONAppr1_EachFreqTime = zeros(dSize,rSize,fSize,tSize); % Capon with summation accross frequency before inversion
%-------------------------------------------------------------------------%


% for nl = 1:lSize
%     for nL = 1:LSize
%         
%         H = (LSearch(nL) - rA + h)*pi/12; % Hour Angle in radians 
%        
%         Y_T = Y(lSearch(nl));
%         Y_T = reshape(Y_T,1,1,6);
%         Z_T = Z(lSearch(nl));
%         Z_T = reshape(Z_T,1,1,6);
%         w = -cosd(dec)*sin(H).*X + sind(dec).*Y_T + cosd(dec)*cos(H).*Z_T;
%         tG = w/c;
%         matchedFilter = exp(-1i*2*pi*fc.*(tG + tI));
%         matchedFilter = permute(matchedFilter,[1 3 2]);
%         
%         matchedFilter = 1/sqrt(6)*squeeze(sum(matchedFilter,2));
% 
%         %%% Matched Filter
%         intensityMF(nl,nL) = sum(sum(ABlT .* matchedFilter));   
% 
%     end
% end
% 
%     mesh(lSearch,LSearch,abs(intensityMF),'FaceColor','interp','LineStyle','none')
% %     xlabel('RA (Hour)')
% %     ylabel('Dec (Degrees)')
%     view(2)
%     colormap('jet')

       
      
              

       
        

      
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Frequency
        %-------------------------------------------------------------------------%
        % fc = linspace(0,100,1024)*1e6;
        f = fc.'; % 50 MHz
        % f = 50*1e6; % 50 MHz
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Signal Model
        %-------------------------------------------------------------------------%
        modelledR = @(l,L) reshape(exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        
        %% First Dervatives
        modelledRdiff_l = @(l,L) reshape((1i*2*pi*dwCygA_l(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        modelledRdiff_L = @(l,L) reshape((1i*2*pi*dwCygA_L(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        
        
        
        
        measuredR = reshape(ABlT.',[],1); 
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Min Function
        %-------------------------------------------------------------------------%
        minFunc = @(l,L) (measuredR - modelledR(l,L)) .* conj(measuredR - modelledR(l,L));
        
        %% First Dervatives
        JacobianMinFunc_l = @(l,L) -modelledRdiff_l(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_l(l,L));
        
        JacobianMinFunc_L = @(l,L) -modelledRdiff_L(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_L(l,L));
        
        JacobianMinFunc = @(l,L) [JacobianMinFunc_l(l,L), JacobianMinFunc_L(l,L)];
        
        %% Second Derivatives
        JacobianMinFunc_l_l = @(l,L) -modelledRdiff_l_l(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            -modelledRdiff_l(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
            -modelledRdiff_l(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_l_l(l,L));
        
        JacobianMinFunc_L_L = @(l,L) -modelledRdiff_L_L(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            -modelledRdiff_L(l,L) .* conj(-modelledRdiff_L(l,L)) + ...
            -modelledRdiff_L(l,L) .* conj(- modelledRdiff_L(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_L_L(l,L));
        
        JacobianMinFunc_l_L = @(l,L) -modelledRdiff_l_L(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            -modelledRdiff_l(l,L) .* conj(-modelledRdiff_L(l,L)) + ...
            -modelledRdiff_L(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_l_L(l,L));
        
        JacobianMinFunc_L_l = @(l,L) -modelledRdiff_L_l(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            -modelledRdiff_L(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
            -modelledRdiff_l(l,L) .* conj(- modelledRdiff_L(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_L_l(l,L));
        
        JacobianMinFunc2 = @(l,L) [JacobianMinFunc_l_l(l,L), JacobianMinFunc_l_L(l,L);...
            JacobianMinFunc_L_l(l,L), JacobianMinFunc_L_L(l,L)];
        %-------------------------------------------------------------------------%
        
        
        
        %% Newton and Secant Search for Two Parameter
        [lL_k,~,all_lL] = newton(minFunc,JacobianMinFunc,[34;-116]);
        % [lL_k,~,all_lL] = newton(JacobianMinFunc,JacobianMinFunc2,[lEst;LEst]);
        
        % [lL_k,~,all_lL] = secant(JacobianMinFunc,[lEst1;LEst1],[lEst2;LEst2]);
        % [lL_k,~,all_lL] = secant(minFunc,[lEst1;LEst1],[lEst2;LEst2]);
        
        % figure;
        % plot(abs(lIn-all_lL(1,:)))
        % figure;
        % plot(all_lL(2,:))
        plot(all_lL)






%% Newton Method
function [x_k,k,all_x] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess: row vector

% x = root of the function

% tol = 10^-4;   % tolerance value
nIter = 1e3; % iteration value

x_k_1 = x_0;    % x(k-1)

all_x = zeros(length(x_0),nIter);
k = 1;
while k <= nIter % norm(x_k - x_k_1) > tol % termination criterion
    
    x_k_1_c = num2cell(x_k_1);
    J = df(x_k_1_c{:});
    f_k_1 = f(x_k_1_c{:});
    
    x_k = x_k_1 - 1*((J' * J) \ J') * f_k_1;
    
    if length(x_0)>1
        all_x(:,k) = x_k;
    else
        all_x(k) = x_k;
    end
    
    x_k_1 = x_k;
    
    k = k + 1;
end

end



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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the GST
%-------------------------------------------------------------------------%
function GST = calculate_GST(hour,minute,second,day,month,year)

% JD2000 = juliandate(2000,1,1,12,0,0);
JD2000 = juliandate(2000,1,1,11,58,55.816);

JD = juliandate(year,month,day);
UT = hour + minute/60 + second/3600;
T = (JD - JD2000)/36525.0;
GST = 6.697374558 + (2400.051336*T)+(0.000025862*T^2)+(UT*1.0027379093);
GST = mod(GST,24);
end
%-------------------------------------------------------------------------%