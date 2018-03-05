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
% Antenna Coordinates
%-------------------------------------------------------------------------%
localLat = 35.247469;
localLong = -116.791509;
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
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LST of the station at the beginning of record
%-------------------------------------------------------------------------%
LST = 100.46 + 0.985647 * numOfDays + localLong + 15*startHour;
LST = mod(LST,360);

% LST of the station during the whole record
% tSize = 40000;
% hRecording = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
% hRecording = (0 : (tSize-1))/3600;

tSize = 18775; % tSize of AB Data
hRecording = (0*12*0.083886 : 12*0.083886 : (tSize-1)*12*0.083886)/3600;

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

thao = haCygnusA - localLong;
thao = mod(thao,360);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a Measurement
%-------------------------------------------------------------------------%
% Antenna Baseline
f = 10e6; % 10 MHz
c = physconst('LightSpeed');
x = 160*f/c;
y = -30*f/c;
z = 3*f/c;

% Fringe Rate
wE = 7.27 * 1e-5;
measuredFR = -wE * (x*cosd(haCygnusA) - y*sind(localLat)*sind(haCygnusA) + z*cosd(localLat)*sind(haCygnusA)) * cosd(decCygnusA);

% Fringe Function
measuredF = -x*cosd(decCygnusA)*sind(haCygnusA) + (z*sind(localLat) + y*cosd(localLat))*sind(decCygnusA) + ...
    (z*cosd(localLat) - y*sind(localLat))*cosd(decCygnusA)*cosd(haCygnusA);

% Correlation
measuredR = cos(2*pi*(-x*cosd(decCygnusA)*sind(haCygnusA) + (z*sind(localLat) + y*cosd(localLat))*sind(decCygnusA) + ...
    (z*cosd(localLat) - y*sind(localLat))*cosd(decCygnusA)*cosd(haCygnusA)));
%-------------------------------------------------------------------------%


SNRdB = 10:5:30; % dB
% SNRdB = 10; % dB
SNR = 10.^(SNRdB/10);

distanceError = zeros(1,length(SNR));

for nS = 1:length(SNR)
    
    noiseVar = var(measuredF)/SNR(nS);
    noise = sqrt(noiseVar)*randn(1,tSize);
    measuredF = measuredF + noise;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solution for Latitude
    %-------------------------------------------------------------------------%
    
    d = decCygnusA;
    
    %% For Fringe Rate
    % H = haCygnusA;
    % [l_newton,k_newton,l_newton_all] = newton(@(l) (measuredFR + wE * (x*cosd(H) - y*sind(l)*sind(H) + z*cosd(l)*sind(H)) * cosd(d)),...
    %     @(l) (wE * cosd(d) * (-y*cosd(l)*sind(H) - z*sind(l)*sind(H))),30);
    
    %% For Fringe Function (latitude)
    % H = haCygnusA;
    % [l_newton,k_newton_l,l_newton_all] = newton(@(l) (measuredF - (-x*cosd(d)*sind(H) + (z*sind(l)+y*cosd(l))*sind(d) + ...
    %     (z*cosd(l)-y*sind(l))*cosd(d)*cosd(H))),...
    %     @(l) -((z*cosd(l)-y*sind(l))*sind(d) + (-z*sind(l)-y*cosd(l))*cosd(d)*cosd(H)),40);
    
    %% For Fringe Function (Longitude)
%     l = localLat;
%     [L_newton,k_newton_L,L_newton_all] = newton(@(L) (measuredF - (-x*cosd(d)*sind(thao+L) + (z*sind(l)+y*cosd(l))*sind(d) + ...
%         (z*cosd(l)-y*sind(l))*cosd(d)*cosd(thao+L))),...
%         @(L) -(-x*cosd(d)*cosd(thao+L) - (z*cosd(l)-y*sind(l))*cosd(d)*sind(thao+L)),30);
    
    %% For Correlation Function (latitude)
    H = haCygnusA;
    [l_newton,k_newton_l,l_newton_all] = newton(@(l) (measuredR - cos(2*pi*(-x*cosd(d)*sind(H) + (z*sind(l)+y*cosd(l))*sind(d) + ...
        (z*cosd(l)-y*sind(l))*cosd(d)*cosd(H)))),...
        @(l) sin(2*pi*(-x*cosd(d)*sind(H) + (z*sind(l)+y*cosd(l))*sind(d) + ...
        (z*cosd(l)-y*sind(l))*cosd(d)*cosd(H))).*(2*pi*((z*cosd(l)-y*sind(l))*sind(d) + (-z*sind(l)-y*cosd(l))*cosd(d)*cosd(H))),...
        10);
    
    %% Secant method
    % [L_secant,k_secant] = secant(@(L) (measuredF + wE * (x*cosd(H) - y*sind(L)*sind(H) + z*cosd(L)*sind(H)) * cosd(decCygnusA)),...
    %     30,48);
    %-------------------------------------------------------------------------%
    
    % For latitude
    distanceError(nS) = distance(localLat,localLong,l_newton,localLong,wgs84Ellipsoid);
    
    % For Longitude
    % distanceError(nS) = distance(localLat,localLong,localLat,L_newton,wgs84Ellipsoid);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Method
%-------------------------------------------------------------------------%
function [x_k,k,x_k_all] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess

% x = root of the function
% k = number of iterations

tol = 10^-7;   % tolerance value
x_k_1 = x_0;
x_k = x_k_1 - (df(x_k_1) * f(x_k_1).') / norm(df(x_k_1))^2;
k = 0;
x_k_all = [];
while abs(x_k - x_k_1) > tol   % termination criterion
    x_k_1 = x_k;
    x_k = x_k_1 - (df(x_k_1) * f(x_k_1).') / norm(df(x_k_1))^2;
    x_k_all = [x_k_all x_k];
    k = k + 1;
end
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Secant Method
%-------------------------------------------------------------------------%
function [x_k,k] = secant(f,x_0,x_1)
% f = nonlinear function
% x_0 = first initial guess
% x_1 = second initial guess
 
% x = root of the function
% k = number of iterations
 
tol = 10^-5;   % tolerance value
x_k_1 = x_0;    % x(k-1)
x_k = x_1;      % x(k)
k = 0;
while abs(f(x_k)) > tol   % termination criterion
    x_k = x_k - f(x_k) * (x_k - x_k_1) / (f(x_k) - f(x_k_1));
    x_k_1 = x_k;
    k = k + 1;
end
end
%-------------------------------------------------------------------------%