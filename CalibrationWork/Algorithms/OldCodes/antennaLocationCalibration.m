fRange = [500:550];
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
% z = -innerNED(3,:); z will be passed as a variable

% plot(x,y,'o') x: towards the east, y: toward the north
% sum((innerNED(:,1)-innerNED(:,2)).^2)^0.5
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baselines and Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = x;
Y = @(z) z*sind(l) + y*cosd(l);
Z = @(z) z*cosd(l) - y*sind(l);

dCygA = decCygnusA;
rCygA = raCygnusA;
HCygA = (LST + h - rCygA).';
wCygA = @(z) -cosd(dCygA)*sin(HCygA*pi/12)*X + sind(dCygA)*Y(z) + cosd(dCygA)*cos(HCygA*pi/12)*Z(z);
dwCygA = sind(dCygA)*sind(l) + cosd(dCygA)*cos(HCygA*pi/12)*cosd(l);

dCasA = decCasA;
rCasA = raCasA;
HCasA = (LST + h - rCasA).';
wCasA = @(z) -cosd(dCasA)*sin(HCasA*pi/12)*X + sind(dCasA)*Y(z) + cosd(dCasA)*cos(HCasA*pi/12)*Z(z);
dwCasA = sind(dCasA)*sind(l) + cosd(dCasA)*cos(HCasA*pi/12)*cosd(l);

%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range and Measured Correlation Data
%-------------------------------------------------------------------------%
fc = (0.09765625*0:0.09765625:0.09765625*1023)*1e6;
c = physconst('LightSpeed');

f = fc(fRange); f = reshape(f,1,1,fSize);
measuredR = ABlT.';
measuredR = measuredR - mean(measuredR);
measuredR = measuredR ./ sqrt(var(measuredR)); % Optional
measuredR = reshape(measuredR,[],1);

%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% radiation Pattern
%-------------------------------------------------------------------------%
% load('LoFASMAntenna60MHz.mat')
% pattern(LoFASMAntenna,60*1e6,0,(H*15).','Type','power');
%-------------------------------------------------------------------------%


% measuredR = real(measuredR)/max(real(measuredR))-mean(real(measuredR));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test measuredR and calculatedR
%-------------------------------------------------------------------------%
% nsF = 0; % 16; 68 ; 135; 69.5; % Best Value is 67
% tI = nsF/(100e6);
% 
% calculatedR = exp(1i*2*pi*f*w/c) * ones(6,1); % *cG_newton;
% plot(real(calculatedR)); % /max(real(calculatedR))))
% figure
% plot(real(measuredR)); %/max(real(measuredR))))

% calculatedR = 1/6 * cos(2*pi*f*(w/c-tI)) * ones(6,1); 
% plot(calculatedR)
% hold on
% plot(measuredR); % /max(real(measuredR)))

% 
% figure; plot(angle(calculatedR))
% hold on
% plot(angle(measuredR))
%-------------------------------------------------------------------------%


%% Minimization on measured R

%% When multiple delays for single frequency
% [T_newton, k_newton, error] = newton(@(T) ((measuredR - exp(1i*2*pi*f*(w/c-T.'))*ones(6,1)) .* ... 
%                                          (conj(measuredR) - exp(-1i*2*pi*f*(w/c-T.'))*ones(6,1))),...
%     @(T) (1i*2*pi*f*exp(1i*2*pi*f*(w/c-T.'))) .* (conj(measuredR) - exp(-1i*2*pi*f*(w/c-T.'))*ones(6,1)) + ...
%     (measuredR - exp(1i*2*pi*f*(w/c-T.'))*ones(6,1)) .* (-1i*2*pi*f*exp(-1i*2*pi*f*(w/c-T.'))), ...
%     zeros(6,1));

%% When single delay for multiple frequency
% Complex
T = 6.7732e-07; % 68.73/(100e6);

gainCygA = repmat(cosd(HCygA*15),fSize,1);
gainCasA = repmat(cosd(HCasA*15),fSize,1);

modelledR = @(z) gainCygA.*reshape(permute(exp(1i*2*pi*f.*(wCygA(z)/c-T)),[1 3 2]),[],6)*ones(6,1) + ...
    gainCasA.*reshape(permute(exp(1i*2*pi*f.*(wCasA(z)/c-T)),[1 3 2]),[],6)*ones(6,1);

modelledRdiff = @(z) gainCygA.*reshape(permute(1i*2*pi/c*f.*dwCygA.*exp(1i*2*pi*f.*(wCygA(z)/c-T)),[1 3 2]),[],6)*ones(6,1) + ...
    gainCasA.*reshape(permute(1i*2*pi/c*f.*dwCasA.*exp(1i*2*pi*f.*(wCasA(z)/c-T)),[1 3 2]),[],6)*ones(6,1);


% Old Version
% minFunc = @(T) (measuredR - modelledR(T)) .* conj(measuredR - modelledR(T));
% JacobianMinFunc = @(T) (-modelledRdiff(T) .* conj(measuredR - modelledR(T))) + ...
%     (measuredR - modelledR(T)) .* conj(-modelledRdiff(T));

% New Version
minFunc = @(z) (measuredR - modelledR(z));
JacobianMinFunc = @(z) -modelledRdiff(z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Jacobian
%-------------------------------------------------------------------------%
z = linspace(-30,30,1000);
mR = zeros(1,length(z));
for n = 1:length(z)
    x_k_1 = z(n);
    J = JacobianMinFunc(x_k_1);
    fm = minFunc(x_k_1);
    % mR(n) = -(J.' * J) \ J.' * fm;
    mR(n) = -imag((J' * J) \ J' * fm);
end
plot(z,mR)
[minR,ind] = min(mR);
zEstimated = z(ind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the minimum
%-------------------------------------------------------------------------%
% [z_newton, k_newton, error] = newton(minFunc,JacobianMinFunc,0);

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
function [x_k,k,error] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess

% x = root of the function
% k = number of iterations

% If root error is wanted to be minimized
tol = 10^2;   % tolerance value
x_k_1 = x_0;
J = df(x_k_1);
x_k = x_k_1 - (J.' * J) \ J.' * f(x_k_1);
k = 0;

% If root error is wanted to be minimized
error = norm(x_k - x_k_1);
while norm(x_k - x_k_1) > tol % (norm(real(x_k - x_k_1)) > tol) || (norm(imag(x_k - x_k_1)) > tol)   % termination criterion
    x_k_1 = x_k;
    J = df(x_k_1);
    % svd(J)
    x_k = x_k_1 - (J.' * J) \ J.' * f(x_k_1);
    
    error = [error norm(x_k - x_k_1)];
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