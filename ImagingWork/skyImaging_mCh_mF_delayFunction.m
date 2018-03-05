% tStart = 1;
% ABlT = ABFiltered(110:820,tStart:end);
[fSize, tSize] = size(ABlTerm);

mainLat = 35 + 14/60 + 50.01/3600;
mainLong = -(116 + 47/60 + 35.78/3600);
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

% mainECEF = lla2ecef([mainLat 360-mainLong  0]);
outLat = 35 + 14/60 + 50.93/3600;
outLong = -(116 + 47/60 + 29.37/3600);
outLatMap = 35.247469;
outLongMap = -116.791509;
% outECEF = lla2ecef([outLat 360-outLong  0]);

% DMS = degrees2dms(angleInDegrees);
% [arclen,az] = distance(mainLatMap,mainLongMap,outLatMap,outLongMap,referenceEllipsoid('earth','m'))

% raCygnusA = (19 + 59/60 + 28.3/3600) * 15;
% decCygnusA = 40 + 44/60 + 02/3600;
% Cygnus A: 19h59m28.3s +40d44m02s
% raCasA = (23 + 23/60 + 27.9/3600) * 15;
% decCasA = 58 + 48/60 + 42/3600;
% Cas A: 23h23m27.9s +58d48m42s

% dec = linspace(-pi/2,pi/2,100);
dec = linspace(0,80,50);
% ra = linspace(-pi/2,pi/2,100);
ra = linspace(-20,50,100);
intensity = zeros(length(dec),length(ra));

recordStart = (5 + 45/60 + 11/3600);
% recordStart = 0;
h = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
h = (h + recordStart) * pi/12;


% fc = (0.09765625*110:0.09765625:0.09765625*820)'*1e6;
% fc = (0.09765625*217:0.09765625:0.09765625*417)'*1e6;
% c = physconst('LightSpeed');

% eDistance = 539.84*cosd(10)*0.3048; % East-West baseline in meters
% nDistance = 539.84*sind(10)*0.3048; % North-South baseline in meters
uDistance = 3; % Zenith baseline in meters

%% New Section
load('LoFASMAntenna60MHz.mat')
% tgApproximate = zeros(length(dec),length(ra),length(h));
tgExact = zeros(length(dec),length(ra),6,length(h));
% svApproximate = zeros(length(dec),2,length(ra),length(fc),length(h));

for nH = 1:length(h)
    for dD = 1:length(dec)
        earthRad    = 6371008.7714; % equatorial radius (meters)
        
        % When Approximate Matched Filter with 1 Main Antenna is Used
%         [xMain,yMain,zMain] = sph2cart(mainLongMap*pi/180+h(nH),mainLatMap*pi/180,earthRad);
        % When Exact Matched Filter with 6 Main Antenna is Used
        [xMain1,yMain1,zMain1] = sph2cart(mainLong1Map*pi/180+h(nH),mainLat1Map*pi/180,earthRad);
        [xMain2,yMain2,zMain2] = sph2cart(mainLong2Map*pi/180+h(nH),mainLat2Map*pi/180,earthRad);
        [xMain3,yMain3,zMain3] = sph2cart(mainLong3Map*pi/180+h(nH),mainLat3Map*pi/180,earthRad);
        [xMain4,yMain4,zMain4] = sph2cart(mainLong4Map*pi/180+h(nH),mainLat4Map*pi/180,earthRad);
        [xMain5,yMain5,zMain5] = sph2cart(mainLong5Map*pi/180+h(nH),mainLat5Map*pi/180,earthRad);
        [xMain6,yMain6,zMain6] = sph2cart(mainLong6Map*pi/180+h(nH),mainLat6Map*pi/180,earthRad);

        [xOut,yOut,zOut] = sph2cart(outLongMap*pi/180+h(nH),outLatMap*pi/180,earthRad+uDistance);
        
        % When Approximate Matched Filter with 1 Main Antenna is Used
%         rxLoFASM = phased.HeterogeneousConformalArray(...
%             'ElementSet',{LoFASMAntenna},...
%             'ElementIndices',[1 2],...
%             'ElementPosition',[xMain xOut; yMain yOut; zMain zOut],...
%             'ElementNormal',[[mainLongMap outLongMap]+h(nH)*180/pi; [mainLatMap outLatMap]-90]); %  [azimuth;elevation]
%         viewArray(rxLoFASM,'ShowIndex','All','ShowNormals',true)
        
        % When Exact Matched Filter with 6 Main Antenna is Used
        rxLoFASM = phased.HeterogeneousConformalArray(...
            'ElementSet',{LoFASMAntenna},...
            'ElementIndices',[1 2 3 4 5 6 7],...
            'ElementPosition',[xMain1 xMain2 xMain3 xMain4 xMain5 xMain6 xOut;...
                               yMain1 yMain2 yMain3 yMain4 yMain5 yMain6 yOut;...
                               zMain1 zMain2 zMain3 zMain4 zMain5 zMain6 zOut],...
            'ElementNormal',[[mainLong1Map mainLong2Map mainLong3Map mainLong4Map mainLong5Map mainLong6Map outLongMap]+h(nH)*180/pi;...
                             [mainLat1Map mainLat2Map mainLat3Map mainLat4Map mainLat5Map mainLat6Map outLatMap]-90]); %  [azimuth;elevation]
%         % viewArray(rxLoFASM,'ShowIndex','All','ShowNormals',true)
        
        % When Delay Parameter is Needed
        delayLoFASM = phased.ElementDelay('SensorArray',rxLoFASM);
        tau = delayLoFASM([ra;dec(dD)*ones(1,length(ra))]);
% %         
% %         % When Approximate Matched Filter with 1 Main Antenna is Used
%         tgApproximate(dD,:,nH) = diff(tau);
%         
%         % When Exact Matched Filter with 6 Main Antenna is Used
        for nA = 1:6
            tgExact(dD,:,nA,nH) = tau(7,:)-tau(nA,:);
        end
        
        % When Steering Vector is Needed
%         svLoFASM = phased.SteeringVector('SensorArray',rxLoFASM);
%         svApproximate(dD,:,:,:,nH) = svLoFASM(fc',[ra;dec(dD)*ones(1,length(ra))]);

        release(delayLoFASM);
%         release(svLoFASM);
    end
end

% nsF = 70;
% thao = nsF/(100e6);
% for dD = 1:length(dec)
%     for rR = 1:length(ra)
% %         tgT = reshape(tgApproximate(dD,rR,:),1,length(h));
%         tgT = reshape(tgExact(dD,rR,:,:),6,length(h));
%         
%         %%% Coherent
% %         phaseComp = exp(-2*pi*1i*fc*(tgT-thao));
%         phaseComp = 0;
%         for nA = 1:6
%             phaseComp = phaseComp + exp(-2*pi*1i*fc*(tgT(nA,:)-thao));
%         end
%         intensity(dD,rR) = sum(sum(ABlT .* phaseComp));
%         
%         
% %         %%% Nonceherent
% %         phaseComp = exp(-2*pi*1i*fc*tgT);
% %         phaseComp = 0;
% %         for nA = 1:6
% %             phaseComp = phaseComp + exp(-2*pi*1i*fc*tgT(nA,:));
% %         end
% %         intensity(dD,rR) = sum(abs(sum(ABlT .* phaseComp , 2)).^2);
%     end
% end

% for r = 1:length(ra)
%     for d = 1:length(dec)
%         tau2 = delayLoFASM([ra(r);dec(d)]);
%         tg = tau2(2)-tau2(1);
%         % tg = (152*sin(h+ra(r))*cos(deg2rad(dec(d))))/c*100e6;
%         % matchedFilter = exp(2*pi*1i*(100:824)'*(tg-75)/1024);
%         matchedFilter = exp(-2*pi*1i*fc*tg);
%         intensity(d,r) = sum(sum(ABtT(:,r) .* matchedFilter));
%     end
% end

% for n = 1:length(rotationAngles)
%     lAxesR = rotz(rotationAngles(n))*azelaxes(0,0); % Positive value in rotz() means rotating counterclocwise
%     eAz = lAxesR(:,2);
%     eEl = lAxesR(:,3);
%     [azAng,~,~] = cart2sph(eAz(1),eAz(2),eAz(3));
%     [~,elAng,~] = cart2sph(eEl(1),eEl(2),eEl(3));
%     azAng = azAng * 180 / pi - 90 + signalAzimuthAngle;
%     elAng = elAng * 180 / pi - 90 + signalElevationAngle;
%     angR = [azAng;elAng];
%     % ang_r = [signalAzimuthAngle;signalElevationAngle];
%     % yR = collectPlaneWave(rxULA,xS,angR,fc);
%     yR = collector(xT,angR);
%     noiseR = 0.00001*(randn(size(yR))+1i*randn(size(yR)));
%     AB1(n) = (yR(:,1) + noiseR(:,1))' * (yR(:,2) + noiseR(:,2)) / length(t);
%     % AB2(:,:,n) = (yR + noiseR)' * (yR + noiseR) / length(t);
% end
%% End of new Section

% for d = 1:length(dec)
%     for r = 1:length(ra)
%         tgEW = ((lxyz(2)*cos(h) - lxyz(1)*sin(h))* sin(ra(r)))/c;
%         tgNS = ((-lxyz(2)*sin(h) - lxyz(1)*cos(h)) * sin(decCenter) * sin(dec(d)))/c;
%         tgUD = (-lxyz(3) * sin(decCenter))/c;
%         tg = tgEW + tgNS + tgUD;
%         phaseComp = exp(-2*pi*1i*fc*(tg-thao));
%         
%         intensity(d,r) = sum(sum(ABtT .* phaseComp));
%     end
% end
% % mesh(ra,dec+decCenter,abs(intensity))

% imagesc(ra,dec,intensity)
% mesh(ra,dec,abs(intensity))
% xlabel('RA (Degrees)')
% ylabel('Dec (Degrees)')
% colormap('jet'); 
% % plot(dec,abs(intensity))