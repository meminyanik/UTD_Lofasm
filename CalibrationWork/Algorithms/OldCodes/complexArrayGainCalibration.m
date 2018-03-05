fRange = [300:550];
% 1 = 0 Hz, 110 = 10 MHz, 217 = 20 MHz, 417 = 40 MHz, 820 = 80 MHz, 1024 = 100 MHz

Whitening = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
% nSAv = 50;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
%-------------------------------------------------------------------------%
if nSAv == 447
    load('ABlTermPZ');
else
    load('AB_Padded.mat')
    ABlTerm = longTermAverage(AB,nSAv);
end

if Whitening
    ABlTerm = whiteningApproach1(ABlTerm);
end

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
uDistance = 3; % Zenith baseline in meters 0.6907 ; % 
outAltMap = distdim(3525,'ft','m') + uDistance; % mainAltMap; % 
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
LST = calculate_LST(outLongMap,5,45,11,21,7,2016);
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

x = innerNED(2,:);
y = innerNED(1,:);
z = -innerNED(3,:); % repmat(0.6907,1,6); %

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
%%% Define Frequency Range and Measured Correlation Data
%-------------------------------------------------------------------------%
fc = (0.09765625*0:0.09765625:0.09765625*1023)*1e6;
c = physconst('LightSpeed');

f = fc(fRange); f = reshape(f,1,1,fSize);
measuredR = ABlT.';
% measuredR = measuredR - mean(measuredR);
% measuredR = measuredR ./ sqrt(var(measuredR)); % Optional
measuredR = reshape(measuredR,[],1);

realMeasuredR = real(measuredR);
imagMeasuredR = imag(measuredR);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Radiation Pattern
%-------------------------------------------------------------------------%
% load('LoFASMAntenna60MHz.mat')

% [pat,azimuth,elevation] = pattern(object,frequency,azimuth,elevation)
% pattern(LoFASMAntenna,60*1e6,0,(HCygA*15).','CoordinateSystem','rectangular','Type','power');
% [pat,az,el] = pattern(LoFASMAntenna,60*1e6,0,0:90,'Type','power');
% [efield,az,el] = pattern(LoFASMAntenna,60*1e6,0,(HCygA*15).','Type','power');
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Symbolics
%-------------------------------------------------------------------------%
% syms m w c f t
% % cG = ones(6,1) * r*exp(1i*t);
% f = abs(m - exp(1i*2*pi*f*(w/c-t)))^2; % *r*exp(1i*t);
% % df_r = diff(f,r);
% df_t = diff(f,t);


%% Minimization on measured R

%% When multiple delays for single frequency
% [T_newton, k_newton, error] = newton(@(T) ((measuredR - exp(1i*2*pi*f*(w/c-T.'))*ones(6,1)) .* ... 
%                                          (conj(measuredR) - exp(-1i*2*pi*f*(w/c-T.'))*ones(6,1))),...
%     @(T) (1i*2*pi*f*exp(1i*2*pi*f*(w/c-T.'))) .* (conj(measuredR) - exp(-1i*2*pi*f*(w/c-T.'))*ones(6,1)) + ...
%     (measuredR - exp(1i*2*pi*f*(w/c-T.'))*ones(6,1)) .* (-1i*2*pi*f*exp(-1i*2*pi*f*(w/c-T.'))), ...
%     zeros(6,1));

%% When single delay for multiple frequency
% Complex

gainCygA = 0.81*repmat(cosd(HCygA*15),1,6,fSize); % 8100
gainCasA = 1.1*repmat(cosd(HCasA*15),1,6,fSize); % 11000

%-------------------------------------------------------------------------%

complexGain = @(T) repmat(1*exp(-1i*2*pi*f.*T),tSize,1,1);
complexGainDiff = @(T) repmat(1*(-1i*2*pi*f).*exp(-1i*2*pi*f.*T),tSize,1,1);

modelledR = @(T) sum(reshape(permute(...
    complexGain(T).*(gainCygA.*exp(1i*2*pi*f.*(wCygA/c)) + gainCasA.*exp(1i*2*pi*f.*(wCasA/c)))...
    ,[1 3 2]),[],6),2);

% Single T
% modelledRdiff = @(T) sum(reshape(permute(...
%     complexGainDiff(T).*(gainCygA.*exp(1i*2*pi*f.*(wCygA/c)) + gainCasA.*exp(1i*2*pi*f.*(wCasA/c)))...
%     ,[1 3 2]),[],6),2);

% Multiple T
modelledRdiff = @(T) reshape(permute(...
    complexGainDiff(T).*(gainCygA.*exp(1i*2*pi*f.*(wCygA/c)) + gainCasA.*exp(1i*2*pi*f.*(wCasA/c)))...
    ,[1 3 2]),[],6);

minFunc = @(T) (measuredR - modelledR(T));
JacobianMinFunc = @(T) -modelledRdiff(T);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% realModelledR = @(T) reshape(permute(...
%     gainCygA.*cos(2*pi*f.*(wCygA/c-T)) + gainCasA.*cos(2*pi*f.*(wCasA/c-T))...
%     ,[1 3 2]),[],6)*ones(6,1);
% realModelledRdiff = @(T) reshape(permute(...
%     gainCygA.*(2*pi*f).*sin(2*pi*f.*(wCygA/c-T)) + gainCasA.*(2*pi*f).*sin(2*pi*f.*(wCasA/c-T))...
%     ,[1 3 2]),[],6)*ones(6,1);
% 
% imagModelledR = @(T) reshape(permute(...
%     gainCygA.*sin(2*pi*f.*(wCygA/c-T)) + gainCasA.*sin(2*pi*f.*(wCasA/c-T))...
%     ,[1 3 2]),[],6)*ones(6,1);
% imagModelledRdiff = @(T) reshape(permute(...
%     gainCygA.*(-2*pi*f).*cos(2*pi*f.*(wCygA/c-T)) + gainCasA.*(-2*pi*f).*cos(2*pi*f.*(wCasA/c-T))...
%     ,[1 3 2]),[],6)*ones(6,1);
% 
% minFunc = @(T) (realMeasuredR - realModelledR(T)).^2 + (imagMeasuredR - imagModelledR(T)).^2;
% JacobianMinFunc = @(T) -2*(realMeasuredR - realModelledR(T)).*realModelledRdiff(T) + ...
%     -2*(imagMeasuredR - imagModelledR(T)).*imagModelledRdiff(T);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the minimum
%-------------------------------------------------------------------------%
% [T_newton, k_newton, allT, error] = newton(minFunc,JacobianMinFunc,10/(100e6)); % Single T
[T_newton, k_newton, allT, error] = newton(minFunc,JacobianMinFunc,repmat(60/(100e6),6,1)); % Multiple T


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test measuredR and modelledR
%-------------------------------------------------------------------------%
% rRtest = modelledR(0);
% plot(real(measuredR(1:504)))
% figure
% plot(real(rRtest(1:504)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Jacobian (Single T)
%-------------------------------------------------------------------------%
% T = 67/(100e6);
% t = linspace(0,128,100)/(100e6);
% mR = zeros(1,length(t));
% for n = 1:length(t)
%     x_k_1 = t(n);
%     J = JacobianMinFunc(x_k_1);
%     fm = minFunc(x_k_1);
%     % mR(n) = -(J.' * J) \ J.' * fm;
%     mR(n) = -imag((J' * J) \ J' * fm);
% end
% plot(t*(100e6),mR)
% [~,ind] = min(mR);
% TEstimated = t(ind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Jacobian (Multiple T)
%-------------------------------------------------------------------------%
% % T = 67/(100e6);
% t = linspace(0,512,1000)/(100e6);
% mR = zeros(1,length(t));
% for n = 1:length(t)
%     x_k_1 = t(n);
%     J = JacobianMinFunc(x_k_1);
%     fm = minFunc(x_k_1);
%     % mR(n) = -(J.' * J) \ J.' * fm;
%     mR(n) = -imag((J' * J) \ J' * fm);
% end
% plot(t*(100e6),mR)
% [~,ind] = min(mR);
% TEstimated = t(ind);


% [T_newton, k_newton, error] = newton(@(T) ((measuredR - reshape(permute(exp(1i*2*pi*f.*(w/c-T)),[1 3 2]),[],6)*ones(6,1)) .* ...
%     (conj(measuredR) - reshape(permute(exp(-1i*2*pi*f.*(w/c-T)),[1 3 2]),[],6)*ones(6,1))),...
%     @(T) (reshape(permute(1i*2*pi*f.*exp(1i*2*pi*f.*(w/c-T)),[1 3 2]),[],6)*ones(6,1)) .* (conj(measuredR) - reshape(permute(exp(-1i*2*pi*f.*(w/c-T)),[1 3 2]),[],6)*ones(6,1)) + ...
%     (measuredR - reshape(permute(exp(1i*2*pi*f.*(w/c-T)),[1 3 2]),[],6)*ones(6,1)) .* (reshape(permute(-1i*2*pi*f.*exp(-1i*2*pi*f.*(w/c-T)),[1 3 2]),[],6)*ones(6,1)), ...
%     0);

% Real
% [T_newton, k_newton1, ~] = newton(@(T) (real(measuredR) - reshape(permute(cos(2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)),...
%     @(T) (reshape(permute(-2*pi*f.*sin(2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)), ...
%     zeros(6,1));
% 
% [T_newton, k_newton2, error] = newton(@(T) (imag(measuredR) - reshape(permute(sin(2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)),...
%     @(T) (reshape(permute(2*pi*f.*cos(2*pi*f.*(w/c-T.')),[1 3 2]),[],6)*ones(6,1)), ...
%     T_newton);
% 
% k_newton = k_newton1 + k_newton2;

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

%% Minimization on real part of measured R
% [T_newton, k_newton, error] = newton(@(T) (real(measuredR) - 1/6*cos(2*pi*f*(w/c-T.'))*ones(6,1)),...
%     @(T) (-2*pi*f*1/6*sin(2*pi*f*(w/c-T.'))),...
%     zeros(6,1));

% [T_newton, k_newton, error] = newton(@(T) (real(measuredR) - (cos(2*pi*f*(w(:,1)/c-T(1))) + ...
%     cos(2*pi*f*(w(:,2)/c-T(2))) + cos(2*pi*f*(w(:,3)/c-T(3))) + ...
%     cos(2*pi*f*(w(:,4)/c-T(4))) + cos(2*pi*f*(w(:,5)/c-T(5))) + cos(2*pi*f*(w(:,6)/c-T(6))))),...
%     @(T) (-2*pi*f * [sin(2*pi*f*(w(:,1)/c-T(1))) , ...
%     sin(2*pi*f*(w(:,2)/c-T(2))) , sin(2*pi*f*(w(:,3)/c-T(3))) , ...
%     sin(2*pi*f*(w(:,4)/c-T(4))) , sin(2*pi*f*(w(:,5)/c-T(5))) , sin(2*pi*f*(w(:,6)/c-T(6)))]),...
%     zeros(6,1));

%% Minimization on imaginary part of measured R
% [T_newton, k_newton, error] = newton(@(T) (imag(measuredR) - sin(2*pi*f*(w/c-T.'))*ones(6,1)),...
%     @(T) (2*pi*f*cos(2*pi*f*(w/c-T.'))),...
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
% Newton Method
%-------------------------------------------------------------------------%
function [x_k,k,allx,error] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess

% x = root of the function
% k = number of iterations

% If root error is wanted to be minimized
% tol = 10^-10; % tolerance value
nIter = 1e4; % iteration value

x_k_1 = x_0;
J = df(x_k_1.');
fm = f(x_k_1.');
% x_k = x_k_1 - (J.' * J) \ J.' * fm;
x_k = x_k_1 - imag((J' * J) \ J' * fm);
k = 1;

% If root error is wanted to be minimized
error = zeros(1,nIter);
allx = zeros(length(x_0),nIter);
while k < nIter % norm(x_k - x_k_1) > tol % termination criterion
    x_k_1 = x_k;
    J = df(x_k_1.');
    fm = f(x_k_1.');
    % x_k = x_k_1 - (J.' * J) \ J.' * fm;
    x_k = x_k_1 - imag((J' * J) \ J' * fm);
    
%     dF = 0.09765625*1e6;
%     delta = 1/dF;
%     x_k = rem(x_k,delta);    
    
    error(k) = norm(f(x_k.'));
    allx(:,k) = x_k;
    k = k + 1;
%     plot(norm(error))
%     drawnow
end

% If function is wanted to be minimized
% error = norm(f(x_k));
% while norm(f(x_k)) > tol   % termination criterion
%     x_k_1 = x_k;
%     J = df(x_k_1);
%     x_k = x_k_1 - (J.' * J) \ J.' * f(x_k_1);
%     error = [error norm(f(x_k))];
%     k = k + 1;
% end
end
%-------------------------------------------------------------------------%


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