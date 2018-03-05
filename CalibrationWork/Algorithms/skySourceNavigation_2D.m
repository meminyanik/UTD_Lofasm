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
localLat = 36.247190;
localLong = 72.793287;
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
tSize = 1e6;
% hRecording = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
hRecording = (0 : (tSize-1))/3600;
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
x = 160;
y = -30;
z = 3;

% Earth Rotation Rate
wE = 7.27 * 1e-5;

measuredF = -wE * (x*cosd(haCygnusA) - y*sind(localLat)*sind(haCygnusA) + z*cosd(localLat)*sind(haCygnusA)) * cosd(decCygnusA);
% measuredF = measuredF.';
%-------------------------------------------------------------------------%

SNRdB = 10:5:30; % dB
SNR = 10.^(SNRdB/10);

distanceError = zeros(1,length(SNR));

for nS = 1:length(SNR)
    
    noiseVar = var(measuredF)/SNR(nS);
    noise = sqrt(noiseVar)*randn(1,tSize);
    
    measuredFN = measuredF + noise;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solution for Latitude
    %-------------------------------------------------------------------------%
    % H = haCygnusA;
    
    [lL_newton, k_newton] = newton(@(lL) (measuredFN + wE * (x*cosd(thao+lL(2)) - y*sind(lL(1))*sind(thao+lL(2)) + z*cosd(lL(1))*sind(thao+lL(2))) * cosd(decCygnusA)),...
        @(lL) (wE * cosd(decCygnusA) * (-y*cosd(lL(1))*sind(thao+lL(2)) - z*sind(lL(1))*sind(thao+lL(2)))),...
        @(lL) (wE * cosd(decCygnusA) * (-x*sind(thao+lL(2)) - y*sind(lL(1))*cosd(thao+lL(2)) + z*cosd(lL(1))*cosd(thao+lL(2)))),...
        [30;60]);
    
    % [L_secant,k_secant] = secant(@(L) (measuredF + wE * (x*cosd(H) - y*sind(L)*sind(H) + z*cosd(L)*sind(H)) * cosd(decCygnusA)),...
    %     10,45);
    %-------------------------------------------------------------------------%
    
    distanceError(nS) = distance(localLat,localLong,lL_newton(1),lL_newton(2),wgs84Ellipsoid);

end

plot(SNRdB,distanceError)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Method
%-------------------------------------------------------------------------%
function [lL_k,k] = newton(f,dfl,dfL,lL_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess

% x = root of the function
% k = number of iterations

% If root error is wanted to be minimized
tol = 10^-10;   % tolerance value
lL_k_1 = lL_0;
J = [dfl(lL_k_1)',dfL(lL_k_1)'];
lL_k = lL_k_1 - (J.' * J) \ J.' * f(lL_k_1)';
k = 0;
while (abs(lL_k(1) - lL_k_1(1)) > tol) || (abs(lL_k(2) - lL_k_1(2)) > tol)   % termination criterion
    lL_k_1 = lL_k;
    J = [dfl(lL_k_1)',dfL(lL_k_1)'];
    lL_k = lL_k_1 - (J.' * J) \ J.' * f(lL_k_1)'; 
    k = k + 1;
end

% If function is wanted to be minimized
% tol = 10^-4;   % tolerance value
% LH = LH_0';
% k = 0;
% while norm(f(LH)) > tol   % termination criterion
%     J = [dfL(LH),dfH(LH)];
%     LH = LH - inv(J.' * J) * J.' * f(LH);
%     k = k + 1;
% end
end
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Secant Method
%-------------------------------------------------------------------------%
% function [L,k] = secant(f,L_0,L_1)
% % f = nonlinear function
% % x_0 = first initial guess
% % x_1 = second initial guess
%  
% % x = root of the function
% % k = number of iterations
%  
% tol = 10^-10;   % tolerance value
% L_k_1 = L_0;    % x(k-1)
% L_k = L_1;      % x(k)
% k = 0;
% while abs(f(L_k)) > tol   % termination criterion
%     L = L_k - f(L_k) * (L_k - L_k_1) / (f(L_k) - f(L_k_1));
%     L_k_1 = L_k;
%     L_k = L;
%     k = k + 1;
% end
% end
%-------------------------------------------------------------------------%