fRange = [500,550];
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

lstOut = mod(outLongMap + GST*15,360);
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
%%% Calculate Total Scanned Hour Duration
%-------------------------------------------------------------------------%
h = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
h = h * 180/12; % in degrees
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baseline from Map Locations
%-------------------------------------------------------------------------%
wgs84 = wgs84Ellipsoid('meters');
% Outtrigger is the reference
% [north1,east1,down1] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
% [north2,east2,down2] = geodetic2ned(mainLat2Map,mainLong2Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
% [north3,east3,down3] = geodetic2ned(mainLat3Map,mainLong3Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
% [north4,east4,down4] = geodetic2ned(mainLat4Map,mainLong4Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
% [north5,east5,down5] = geodetic2ned(mainLat5Map,mainLong5Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
% [north6,east6,down6] = geodetic2ned(mainLat6Map,mainLong6Map,mainAltMap,outLatMap,outLongMap,outAltMap,wgs84);
% 
% x = [east1, east2, east3, east4, east5, east6];
% y = [north1, north2, north3, north4, north5, north6];
% z = -1*[down1, down2, down3, down4, down5, down6];

%-------------------------------------------------------------------------%
% Test Lengths at Inner Ring
% [north12,east12,down12] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,mainLat2Map,mainLong2Map,mainAltMap,wgs84);
% [north13,east13,down13] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,mainLat3Map,mainLong3Map,mainAltMap,wgs84);
% [north14,east14,down14] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,mainLat4Map,mainLong4Map,mainAltMap,wgs84);
% [north15,east15,down15] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,mainLat5Map,mainLong5Map,mainAltMap,wgs84);
% [north16,east16,down16] = geodetic2ned(mainLat1Map,mainLong1Map,mainAltMap,mainLat6Map,mainLong6Map,mainAltMap,wgs84);
% sqrt(north12^2+east12^2)
% sqrt(north13^2+east13^2)
% sqrt(north14^2+east14^2)
% sqrt(north15^2+east15^2)
% sqrt(north16^2+east16^2)
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Baseline from Map Locations
%-------------------------------------------------------------------------%
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

% sqrt(sum(innerNED(:,1).^2))
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Theoric Delay
%-------------------------------------------------------------------------%
l = outLatMap;
L = outLongMap;

X = x;
Y = z*sind(l) + y*cosd(l);
Z = z*cosd(l) - y*sind(l);

d = decCygnusA;
r = raCygnusA * 15;
H = (lstOut + h - r).';

w = -cosd(d)*sind(H)*X + sind(d)*Y + cosd(d)*cosd(H)*Z;
c = physconst('LightSpeed');

fc = (0.09765625*0:0.09765625:0.09765625*1023)*1e6;
fRange = [500,550];
f = fc(fRange); f = reshape(f,1,1,length(f));
measuredR = ABlT(fRange,:)'; measuredR = reshape(measuredR,[],1);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test measuredR and calculatedR
%-------------------------------------------------------------------------%
% nsF = 0; % 16; 68 ; 135; 69.5; % Best Value is 67
% tI = nsF/(100e6);
% 
calculatedR = exp(1i*2*pi*f*w/c) * ones(6,1); % *cG_newton;
plot(real(calculatedR/max(real(calculatedR))))
hold on
plot(real(measuredR/max(real(measuredR))))
% 
% figure; plot(angle(calculatedR))
% hold on
% plot(angle(measuredR))
%-------------------------------------------------------------------------%

% Calibrate for Instrument Delay
nsF = 67; % 16; 68 ; 135; 69.5; % Best Value is 67
tI = nsF/(100e6);
% cG_0 = ones(6,1)*exp(-1i*2*pi*f*tI);

% syms m w c f t
% % cG = ones(6,1) * r*exp(1i*t);
% f = m - exp(1i*2*pi*f*(w/c-t)); % *r*exp(1i*t);
% % df_r = diff(f,r);
% df_t = diff(f,t);

% [T_newton, k_newton, error] = newton(@(T) (measuredR - exp(1i*2*pi*f*(w/c-T.'))*ones(6,1)),...
%     @(T) (1i*2*pi*f*exp(1i*2*pi*f*(w/c-T.'))),...
%     zeros(6,1));

[T_newton, k_newton, error] = newton(@(T) (measuredR - (exp(1i*2*pi*f*(w(:,1)/c-T(1))) + ...
    exp(1i*2*pi*f*(w(:,2)/c-T(2))) + exp(1i*2*pi*f*(w(:,3)/c-T(3))) + ...
    exp(1i*2*pi*f*(w(:,4)/c-T(4))) + exp(1i*2*pi*f*(w(:,5)/c-T(5))) + exp(1i*2*pi*f*(w(:,6)/c-T(6))))),...
    @(T) (1i*2*pi*f * [exp(1i*2*pi*f*(w(:,1)/c-T(1))) , ...
    exp(1i*2*pi*f*(w(:,2)/c-T(2))) , exp(1i*2*pi*f*(w(:,3)/c-T(3))) , ...
    exp(1i*2*pi*f*(w(:,4)/c-T(4))) , exp(1i*2*pi*f*(w(:,5)/c-T(5))) , exp(1i*2*pi*f*(w(:,6)/c-T(6)))]),...
    zeros(6,1));


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
tol = 10^-8;   % tolerance value
x_k_1 = x_0;
J = df(x_k_1);
x_k = x_k_1 - (J.' * J) \ J.' * f(x_k_1);
k = 0;

% If root error is wanted to be minimized
% error = abs(x_k - x_k_1);
% while abs(x_k - x_k_1) > tol % (norm(real(x_k - x_k_1)) > tol) || (norm(imag(x_k - x_k_1)) > tol)   % termination criterion
%     x_k_1 = x_k;
%     J = df(x_k_1);
%     x_k = x_k_1 - (J.' * J) \ J.' * f(x_k_1);
%     error = [error abs(x_k - x_k_1)];
%     k = k + 1;
%     % plot(norm(error))
%     % drawnow
% end

% If function is wanted to be minimized
error = norm(f(x_k));
while norm(f(x_k)) > tol   % termination criterion
    x_k_1 = x_k;
    J = df(x_k_1);
    x_k = x_k_1 - (J.' * J) \ J.' * f(x_k_1);
    error = [error norm(f(x_k))];
    k = k + 1;
end
end
%-------------------------------------------------------------------------%