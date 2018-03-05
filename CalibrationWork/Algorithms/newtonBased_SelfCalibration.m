% --------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
currentPath = pwd();
addpath(genpath([currentPath,'/../../../Algorithms']));
addpath(genpath([currentPath,'/../../../RecordedData']));

%% Initial Variables
fRange = 300:700;
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

SimulationMode = false;
Whitening = true;
antennaGainFactor = true;

sourceModelMode = 3; % 1: CygA, 2: Cas A, 3: CygA+CasA

stationID = 1;
arrayConfig = 1;
% 0: inner ring - outtrigger
% 1: outer ring - outtrigger
% 2: inner ring - outer ring
%-------------------------------------------------------------------------%

%% Data Record Time
% Data Record Start Time [h,m,s,d,m,y]
if (stationID == 1)
    dataRecordTime = [20,30,20,3,11,2017]; % 16.89 LST
elseif (stationID == 4)
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
RlT = reshape(RlT,1,[]);
%-------------------------------------------------------------------------%

%% Get Antenna Coordinates
getAntennaCoordinates;
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

% rotNed2BodyMatlab = angle2dcm( (innerRingRotation)*pi/180, 0, 0 );

innerNED = [4.41/2*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
    4.41/2*[-1, -2, -1, 1, 2, 1];...
    zeros(1,6)];

outerNED = [7.6383/2*sqrt(3)*[-1, 0, 1, 1, 0, -1];...
    7.6383/2*[-1, -2, -1, 1, 2, 1];...
    zeros(1,6)];

% plot(innerNED(2,:),innerNED(1,:),'.')

%% When Google Map Variables are used
% innerRingRotation = 90-azimuth(inMainLat2Map,inMainLong2Map,inMainLat5Map,inMainLong5Map,wgs84); % Counter-clockwise - Degrees
% % outerRingRotation = 90-azimuth(outMainLat2Map,outMainLong2Map,outMainLat5Map,outMainLong5Map,wgs84); % Counter-clockwise - Degrees
% outerRingRotation = innerRingRotation + 30;
% rotInNed2Body = [cosd(innerRingRotation) sind(innerRingRotation) 0; ...
%     -sind(innerRingRotation) cosd(innerRingRotation) 0; ...
%     0 0 1];
% rotOutNed2Body = [cosd(outerRingRotation) sind(outerRingRotation) 0; ...
%     -sind(outerRingRotation) cosd(outerRingRotation) 0; ...
%     0 0 1];

%% When Ring rotation is a variable
thIn_Est = 90-azimuth(inMainLat2Map,inMainLong2Map,inMainLat5Map,inMainLong5Map,wgs84); % Counter-clockwise - Degrees
if (arrayConfig == 0)
    x = @(dx,th) [-sind(th) cosd(th) 0]*innerNED + eastMO + dx;
    y = @(dy,th) [cosd(th) sind(th) 0] * innerNED + northMO + dy;
    z = @(dz,th) -[0 0 1] * innerNED + downMO + dz;
    th_Est = thIn_Est;
elseif (arrayConfig == 1)
    x = @(dx,th) [-sind(th) cosd(th) 0] * outerNED + eastMO + dx;
    y = @(dy,th) [cosd(th) sind(th) 0] * outerNED + northMO + dy;
    z = @(dz,th) -[0 0 1] * outerNED + downMO + dz;
    th_Est = thIn_Est + 30;
end

% plot(x(0,th_Est),y(0,th_Est),'o')
% x: towards the east, y: toward the north
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baselines and Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = @(dx,th) x(dx,th);
Y = @(dy,dz,th) z(dz,th)*sind(l) + y(dy,th)*cosd(l);
Z = @(dy,dz,th) z(dz,th)*cosd(l) - y(dy,th)*sind(l);

dCygA = decCygnusA;
rCygA = raCygnusA;
HCygA = (LST + h - rCygA).';
wCygA = @(dx,dy,dz,th) -cosd(dCygA)*sin(HCygA*pi/12)*X(dx,th) + sind(dCygA)*Y(dy,dz,th) + cosd(dCygA)*cos(HCygA*pi/12)*Z(dy,dz,th);

dCasA = decCasA;
rCasA = raCasA;
HCasA = (LST + h - rCasA).';
wCasA = @(dx,dy,dz,th) -cosd(dCasA)*sin(HCasA*pi/12)*X(dx,th) + sind(dCasA)*Y(dy,dz,th) + cosd(dCasA)*cos(HCasA*pi/12)*Z(dy,dz,th);
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
if antennaGainFactor
    gainCygA = cosd(decCygnusA-outLatMap)*0.81*cosd(HCygA.'*15); % 8100
    gainCasA = cosd(decCasA-outLatMap)*1.1*cosd(HCasA.'*15); % 11000
    gainCygA(gainCygA<=0) = 0;
    gainCasA(gainCasA<=0) = 0;
    
    gainCygA = repmat(gainCygA,fSize,1);
    gainCasA = repmat(gainCasA,fSize,1);
else
    gainCygA = 1;
    gainCasA = 1;
end

tGcygA = @(dx,dy,dz,th) reshape(wCygA(dx,dy,dz,th),1,tSize,6)/c;
matchedFilterCygA = @(dx,dy,dz,th,tI) reshape(gainCygA .* sum(exp(-1i*2*pi*fc.*(tGcygA(dx,dy,dz,th)+tI)),3),[],1);
matchedFilterCygA_diff_tI = @(dx,dy,dz,th,tI) reshape(gainCygA .* ...
    sum((-1i*2*pi*fc).*exp(-1i*2*pi*fc.*(tGcygA(dx,dy,dz,th)+tI)),3),[],1);
% matchedFilterCygA = @(th,tI) gainCygA .* sum(exp(-1i*2*pi*fc.*(tGcygA(0,0,0,th)+tI)),3);

tGcasA = @(dx,dy,dz,th) reshape(wCasA(dx,dy,dz,th),1,tSize,6)/c;
matchedFilterCasA = @(dx,dy,dz,th,tI) reshape(gainCasA .* sum(exp(-1i*2*pi*fc.*(tGcasA(dx,dy,dz,th)+tI)),3),[],1);
matchedFilterCasA_diff_tI = @(dx,dy,dz,th,tI) reshape(gainCasA .* ...
    sum((-1i*2*pi*fc).*exp(-1i*2*pi*fc.*(tGcasA(dx,dy,dz,th)+tI)),3),[],1);
% matchedFilterCasA = @(th,tI) gainCasA .* sum(exp(-1i*2*pi*fc.*(tGcasA(0,0,0,th)+tI)),3);
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Peak Intensity for a Specific Location
%---------------------------------------------------------------------%
cygAPeakIntensity = @(dx,dy,dz,th,tI) RlT*matchedFilterCygA(dx,dy,dz,th,tI);
% cygAPeakIntensity = @(th,tI) sum(sum(RlT.*matchedFilterCygA(th,tI)));
% figure; plot(fc,unwrap(angle(cygAPeakIntensity)));
casAPeakIntensity = @(dx,dy,dz,th,tI) RlT*matchedFilterCasA(dx,dy,dz,th,tI);
% casAPeakIntensity = @(th,tI) sum(sum(RlT.*matchedFilterCasA(th,tI)));
% figure; plot(fc,unwrap(angle(casAPeakIntensity)));
cygAcasAPeakIntensity = @(dx,dy,dz,th,tI) RlT*(matchedFilterCygA(dx,dy,dz,th,tI)+matchedFilterCasA(dx,dy,dz,th,tI));

switch sourceModelMode
    case 1
        peakIntensity = @(dx,dy,dz,th,tI) cygAPeakIntensity(dx,dy,dz,th,tI);
    case 2
        peakIntensity = @(dx,dy,dz,th,tI) casAPeakIntensity(dx,dy,dz,th,tI);
    case 3
        peakIntensity = @(dx,dy,dz,th,tI) cygAcasAPeakIntensity(dx,dy,dz,th,tI);
end

%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scan Peak For Instrument Delay
%---------------------------------------------------------------------%
dEst = linspace(100,300,100); % In Meters
tIScan = -dEst/c;
peakScan = zeros(1,length(tIScan));
for n = 1:length(tIScan)
    peakScan(n) = peakIntensity(0,0,0,th_Est,tIScan(n));
end
figure; plot(tIScan,abs(peakScan))
[~,ind] = max(abs(peakScan));
tI_Estimated = tIScan(ind);

% tI_Estimated = -5.927576089874496e-07; % For L1_AC Data
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scan Peak For Z-Offset Rotation
%---------------------------------------------------------------------%
dzScan = linspace(-10,10,100); % In Meters
peakScan = zeros(1,length(dzScan));
for n = 1:length(dzScan)
    peakScan(n) = peakIntensity(0,0,dzScan(n),th_Est,tI_Estimated);
end
figure; plot(dzScan,abs(peakScan))
[~,ind] = max(abs(peakScan));
dz_Estimated = dzScan(ind);

% dz_Estimated = -0.505050505050505; % For L1_AC Data
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scan Peak For Y-Offset Rotation
%---------------------------------------------------------------------%
% dyScan = linspace(-10,10,100); % In Meters
% peakCyg = zeros(1,length(dyScan));
% for n = 1:length(dyScan)
%     peakCyg(n) = cygAPeakIntensity(x_Estimated,dyScan(n),0,th_Est,tI_Estimated);
% end
% figure; plot(dyScan,abs(peakCyg))
% [~,ind] = max(abs(peakCyg));
% dy_Estimated = dyScan(ind);
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scan Peak For X-Offset Rotation
%---------------------------------------------------------------------%
% dxScan = linspace(-10,10,100); % In Meters
% peakCyg = zeros(1,length(dxScan));
% for n = 1:length(dxScan)
%     peakCyg(n) = cygAPeakIntensity(dxScan(n),0,0,th_Est,tI_Estimated);
% end
% figure; plot(dxScan,abs(peakCyg))
% [~,ind] = max(abs(peakCyg));
% dx_Estimated = dxScan(ind);
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scan Peak For Ring Rotation
%---------------------------------------------------------------------%
% thScan = linspace(th_Est-5,th_Est+10,100); % In Meters
% peakCyg = zeros(1,length(thScan));
% for n = 1:length(thScan)
%     peakCyg(n) = cygAPeakIntensity(x_Estimated,0,z_Estimated,thScan(n),tI_Estimated);
% end
% plot(thScan,abs(peakCyg))
% [~,ind] = max(abs(peakCyg));
% th_Estimated = thScan(ind);
%---------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Min Function
%-------------------------------------------------------------------------%
measuredR = RlT.';

switch sourceModelMode
    case 1
        modelledR = @(tI) conj(matchedFilterCygA(0,0,dz_Estimated,th_Est,tI));
        modelledR_diff_tI = @(tI) conj(matchedFilterCygA_diff_tI(0,0,dz_Estimated,th_Est,tI));
    case 2
        modelledR = @(tI) conj(matchedFilterCasA(0,0,dz_Estimated,th_Est,tI));
        modelledR_diff_tI = @(tI) conj(matchedFilterCasA_diff_tI(0,0,dz_Estimated,th_Est,tI));
    case 3
        modelledR = @(tI) conj(matchedFilterCygA(0,0,dz_Estimated,th_Est,tI)) + ...
            conj(matchedFilterCasA(0,0,dz_Estimated,th_Est,tI));
        modelledR_diff_tI = @(tI) conj(matchedFilterCygA_diff_tI(0,0,dz_Estimated,th_Est,tI)) + ...
            conj(matchedFilterCasA_diff_tI(0,0,dz_Estimated,th_Est,tI));
end

minFunc = @(tI) (measuredR - modelledR(tI)) .* conj(measuredR - modelledR(tI));

%% Define Jacobian
JacobianMinFunc = @(tI) -modelledR_diff_tI(tI) .* conj(measuredR - modelledR(tI)) + ...
    (measuredR - modelledR(tI)) .* conj(-modelledR_diff_tI(tI));


%% Find the minimum
%-------------------------------------------------------------------------%
[tI_k,~,all_tI] = newton(minFunc,JacobianMinFunc,tI_Estimated);
tI_Estimated = tI_k;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Method
%-------------------------------------------------------------------------%
function [x_k,k,all_x] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess: row vector

% x = root of the function
% k = number of iterations
% all_x = all roots at each iteraion

% tol = 10^-12;   % tolerance value (when tol is used termination criterion)
nIter = 1e2; % iteration value (when nIter is used termination criterion)

x_k_1 = x_0;    % x(k-1)

all_x = zeros(length(x_0),nIter);
k = 1;
while k <= nIter % norm(x_k - x_k_1) > tol
    
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


%% When multiple delays for multiple frequency
% [T_newton, k_newton, error] = newton(@(T) ((measuredR - reshape(permute(exp(1i*2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)) .* ...
%     (conj(measuredR) - reshape(permute(exp(-1i*2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1))),...
%     @(T) (reshape(permute(1i*2*pi*f.*exp(1i*2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)) .* (conj(measuredR) - reshape(permute(exp(-1i*2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)) + ...
%     (measuredR - reshape(permute(exp(1i*2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)) .* (reshape(permute(-1i*2*pi*f.*exp(-1i*2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)), ...
%     zeros(6,1));

% [T_newton, k_newton, error] = newton(@(T) (measuredR - (exp(1i*2*pi*f*(w(:,1)/c-T(1))) + ...
%     exp(1i*2*pi*f*(w(:,2)/c-T(2))) + exp(1i*2*pi*f*(w(:,3)/c-T(3))) + ...
%     exp(1i*2*pi*f*(w(:,4)/c-T(4))) + exp(1i*2*pi*f*(w(:,5)/c-T(5))) + exp(1i*2*pi*f*(w(:,6)/c-T(6))))),...
%     @(T) (1i*2*pi*f * [exp(1i*2*pi*f*(w(:,1)/c-T(1))) , ...
%     exp(1i*2*pi*f*(w(:,2)/c-T(2))) , exp(1i*2*pi*f*(w(:,3)/c-T(3))) , ...
%     exp(1i*2*pi*f*(w(:,4)/c-T(4))) , exp(1i*2*pi*f*(w(:,5)/c-T(5))) , exp(1i*2*pi*f*(w(:,6)/c-T(6)))]),...
%     zeros(6,1));

%% minimization with fmincon
% measuredR = measuredR(1);
% w = w(1,:);
% fun = @(T)(real(measuredR) - 1/6*cos(2*pi*f*(w/c-T))*ones(6,1));
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = zeros(1,6);
% ub = 1e-4 * ones(1,6);
% T0 = zeros(1,6);
% T_fmincon = fmincon(fun,T0,A,b,Aeq,beq,lb,ub);
% optimtool

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Further Analysis:
% tI = T_newton; % 0.0068; 0.0011
% dF = 0.09765625*1e6;
% delta = 1/dF;
% tI = rem(tI,delta);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Symbolics
%-------------------------------------------------------------------------%
% syms m w c f t
% % cG = ones(6,1) * r*exp(1i*t);
% f = abs(m - exp(1i*2*pi*f*(w/c-t)))^2; % *r*exp(1i*t);
% % df_r = diff(f,r);
% df_t = diff(f,t);