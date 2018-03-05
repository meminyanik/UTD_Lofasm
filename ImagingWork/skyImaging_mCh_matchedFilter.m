tStart = 1;

% RFiltered = RFiltered / 1e6;

% AAlT = reshape(RFiltered(1,1,:,:),1024,504);
% ABlT = reshape(RFiltered(2,1,:,:),1024,504);
% BBlT = reshape(RFiltered(2,2,:,:),1024,504);

% AAlT = BBlTerm;
ABlT = ABlTerm / 1024;
% BBlT = AAlTerm;

% AAlT = AAlT(110:820,tStart:end);
ABlT = ABlT(110:820,tStart:end);
% BBlT = BBlT(110:820,tStart:end);

% AAlT = AAlT(217:417,tStart:end);
% ABlT = ABlT(217:417,tStart:end);
% BBlT = BBlT(217:417,tStart:end);

[fSize, tSize] = size(ABlT);

mainLat = 35 + 14/60 + 50.01/3600;
mainLong = -(116 + 47/60 + 35.78/3600);
mainLatMap = 35.247227;
mainLongMap = -116.793272;
% mainECEF = lla2ecef([mainLat 360-mainLong  0]);

outLat = 35 + 14/60 + 50.93/3600;
outLong = -(116 + 47/60 + 29.37/3600);
outLatMap = 35.247469;
outLongMap = -116.791509;
% outECEF = lla2ecef([outLat 360-outLong  0]);

raCygnusA = (19 + 59/60 + 28.3/3600) * 15;
decCygnusA = 40 + 44/60 + 02/3600;
% Cygnus A: 19h59m28.3s +40d44m02s
raCasA = (23 + 23/60 + 27.9/3600) * 15;
decCasA = 58 + 48/60 + 42/3600;
% Cas A: 23h23m27.9s +58d48m42s

% dec = linspace(-pi/2,pi/2,100);
dec = linspace(20,70,50);
% ra = linspace(-pi/2,pi/2,100);
ra = linspace(-80,80,100); % or -80 to -10
intensity = zeros(length(dec),length(ra));


recordStart = (5 + 45/60 + 11/3600);
% recordStart = 0;
h = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
h = (h + recordStart) * pi/12;

fc = (0.09765625*110:0.09765625:0.09765625*820)'*1e6;
% fc = (0.09765625*110:0.09765625:0.09765625*820)'/100;
c = physconst('LightSpeed');

[arclen,az] = distance(mainLatMap,mainLongMap,outLatMap,outLongMap,referenceEllipsoid('earth','m')); % same with 'wgs84'
eDistance = arclen*cosd(90-az);
nDistance = arclen*sind(90-az);
% eDistance = 539.84*cosd(10)*0.3048; % East-West baseline in meters
% nDistance = 539.84*sind(10)*0.3048; % North-South baseline in meters
uDistance = 3; % Zenith baseline in meters

X = eDistance*cosd(h*180/pi+mainLongMap);
Y = nDistance;
Z = eDistance*sind(h*180/pi+mainLongMap);

ABlT = ABlT .* exp(2*pi*1i*fc/c*Z);

% xyz = zeros(3,length(h));
% for nL = 1:length(h)
%     mainLongT = mainLongMap+h(nL)*180/pi;
%     outLongT = outLongMap+h(nL)*180/pi;
%     
%     enu = [eDistance;nDistance;uDistance];
%     rotM1 = [0 -sin(mainLatMap) cos(mainLatMap);1 0 0;0 cos(mainLatMap) sin(mainLatMap)];
%     xyz(:,nL) = rotM1 * enu;
% end
% X = eDistance;
% Y = uDistance*sind(mainLatMap) + nDistance*cosd(mainLatMap);
% Z = uDistance*cosd(mainLatMap) - nDistance*sind(mainLatMap);

nsF = 70;
thao = nsF/(100e6);

for d = 1:length(dec)
    for r = 1:length(ra)
        % hT = h + ra(r)*pi/180;
        % tgEW = ((lxyz(2)*cos(h) - lxyz(1)*sin(h))* sin(ra(r)))/c;
        % tgNS = ((-lxyz(2)*sin(h) - lxyz(1)*cos(h)) * sin(decCenter) * sin(dec(d)))/c;
        % tgUD = (-lxyz(3) * sin(decCenter))/c;
        % tg = tgEW + tgNS + tgUD;
        % phaseComp = exp(-2*pi*1i*fc*(tg-thao));
%         tgT = reshape(tg(d,r,:),1,length(h));
%         % phaseComp = exp(-2*pi*1i*fc*(tgT-thao));
%         phaseComp = exp(-2*pi*1i*fc*tgT);
%         intensity(d,r) = sum(abs(sum(ABtT .* phaseComp , 2)).^2);
        % tg = (-X*cosd(dec(d))*sin(hT) + Y*sind(dec(d)) + Z*cosd(dec(d))*cos(hT))/c;
        tg = (X*sind(dec(d))*cosd(ra(r)) + Y*sind(dec(d))*sin(ra(r)))/c;
        phaseComp = exp(-2*pi*1i*fc*(tg-thao)); % 2*pi*fc/c is k (wavelength)
        intensity(d,r) = sum(sum(ABlT .* phaseComp));
        % intensity(d,r) = sum(abs(sum(ABtT .* phaseComp , 2)).^2);
    end
end
% mesh(ra,dec+decCenter,abs(intensity))
% mesh(ra,dec,abs(intensity))
mesh(ra,dec,abs(intensity))