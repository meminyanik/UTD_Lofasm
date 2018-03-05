tic
% % % load('AB_054511_084543.mat')
% load('AB_054511_105606.mat')
t1 = toc;

ABt = ABlTerm(110:820,:);
[fSize, tSize] = size(ABt);

% dec = linspace(0.3,1.2,60);
% decCenter = 34.200833*pi/180; % LoFASM IV
decCenter = 90*pi/180;

decCassopedia = 58.8000; % in degrees
decCygnus = 40.7339; % in degrees

dec = linspace(-pi/2,pi/2,100);
% ra = (1.5:-0.05:-3)*pi/12;
% dec = 10:1:60; % declination in degrees
% ra = linspace(-1.5*pi/12,3*pi/12,100);
% ra=(1.5:-0.05:-3)*pi/12;
ra = linspace(-pi/2,pi/2,100);
intensity = zeros(length(dec),length(ra),fSize);

% h = (0.083886*100000:0.083886:0.083886*100020)/3600; % 1 Hour
% h = (0.083886*0:0.083886:0.083886*221723)/3600;
% h = 0.083886*42914/3600;
% h = (10*0:10:10*1863)/3600;
h = (0*447*0.083886 : 447*0.083886 : 495*447*0.083886)/3600;
h = h - 1.25;
h = h * pi/12;

% h = linspace(-1.5, 5-1.5+1/6, tSize+1)*pi/12;
% h = h(1:end-1);

% f index for 10 MHz = 110
% f index for 80 MHz = 820
fc = (0.09765625*110:0.09765625:0.09765625*820)'*1e6;
% fc = flip(fc);
% fc = (0.09765625*600)*1e6;
c = physconst('LightSpeed');

% lat = 34.200833*pi/180;
% enu = [0 -152;0 0;0 0];
% rotM1 = [0 -sin(lat) cos(lat);1 0 0;0 cos(lat) sin(lat)];
% xyz = rotM1 * enu;
% lxyz = xyz(:,1) - xyz(:,2);
% t2 = toc;

nsF = 70;
%%%% thao = 1/(0.09765625*nsF*1e6);
thao = nsF/(100e6);
        
for d = 1:length(dec)
    for r = 1:length(ra)
        % hT = h + ra(r);
        % tg = (152*(cos(hT) + sin(dec(d))*sin(hT)))/c;
        % tg = (152 * sin(hT)*cos(dec(d)))/c;
        
%         tg = (152*sin(h+ra(r))*cos(deg2rad(dec(d))))/c*100e6;
%         phaseComp = exp(2*pi*1i*[100:824]'*(tg-75)/1024);
        
        
                tgEW = (152 * cos(h) * sin(ra(r)))/c;
                tgNS = (152 * sin(h) * sin(decCenter) * sin(dec(d)))/c;
                tg = tgEW + tgNS;
                phaseComp = exp(-2*pi*1i*fc*(tg-thao));
        % tg = (152*(cos(hT) + sin(dec(d))*sin(hT) - sin(hT)*cos(dec(d))))/c;
        
        
        %         hT = reshape(hT,1,1,length(hT));
        %         zeroT = zeros(1,1,length(hT));
        %         decT = dec(d) * ones(1,1,length(hT));
        %         rotM2 = [sin(hT) cos(hT) zeroT;...
        %             -sin(dec(d))*cos(hT) sin(dec(d))*sin(hT) cos(decT);...
        %             cos(dec(d))*cos(hT) -cos(dec(d))*sin(hT) sin(decT)];
        %         uvw = sum(rotM2 .* lxyz',2);
        %         tg = sum(uvw) / c;
        %         tg = reshape(tg,1,length(hT));
        %         % tg = uvw(3,:) / c;
        % phaseComp = exp(2*pi*1i*fc*(tg + thao));
        
        % tg = (152*sin(ht)*cos(dec(d)))/c*100e6;
        % phaseComp = exp(2*pi*1i*(110:820)'*(tg-75)/1024);
        
        intensity(d,r,:) = sum(ABt .* phaseComp,2);
        % intensity(d,r) = ABt * exp(2*pi*1i*fc*tg);
        
        %         hT = reshape(h,1,1,length(h));
        %         zeroT = zeros(1,1,length(h));
        %         decT = dec(d) * ones(1,1,length(h));
        %         rotM2 = [sin(hT) cos(hT) zeroT;...
        %             -sin(dec(d))*cos(hT) sin(dec(d))*sin(hT) cos(decT);...
        %             cos(dec(d))*cos(hT) -cos(dec(d))*sin(hT) sin(decT)];
        
        %         uvw = sum(rotM2 .* lxyz',2) .* (fc/c);
        %         AB2 = reshape(ABt,[1 length(fc)*length(h)]);
        %         uArray = reshape(uvw(1,:,:),[1 length(fc)*length(h)]);
        %         vArray = reshape(uvw(2,:,:),[1 length(fc)*length(h)]);
        %         wArray = reshape(uvw(3,:,:),[1 length(fc)*length(h)]);
        %         % AB2 = AB2 .* exp(-2*pi*1i*uArray) .* exp(-2*pi*1i*vArray) .* exp(-2*pi*1i*wArray);
        %         n = sqrt(1 - (cos(ra(r))^2 + cos(dec(d))^2));
        %         AB2 = AB2 .* exp(-2*pi*1i*cos(ra(r))*uArray) .* exp(-2*pi*1i*cos(dec(d))*vArray) .* exp(-2*pi*1i*n*wArray);
        %         % AB2 = AB2 .* exp(-2*pi*1i*sin(ra(r))*uArray) .* exp(-2*pi*1i*sin(dec(d))*vArray);
        
        %         intensity(d,r) = sum(AB2);
        
        
        %         mesh(ra,dec,abs(sum(intensity,3)))
        %         hold on;
        %         pause(0.01);
        t3 = toc;
    end
end
mesh(ra,dec+decCenter,abs(sum(intensity,3)))
% mesh(ra,dec,fliplr(sum(abs(intensity),3)))
t4 = toc;
