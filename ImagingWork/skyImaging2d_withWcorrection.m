% tic
% % % load('AB_054511_084543.mat')
% load('AB_054511_105606.mat')
% % % t1 = toc;
AB1 = AB(110:820,:);
% % AB = AB(600:650,:);
% % AB = AB(500,:);
% % AB = AB(110:820,:);
% [fSize, tSize] = size(AB);
% % % AB1 = zeros(tSize, 2*fSize-1);
% AB1 = zeros(fSize, tSize);
% for n = 1:tSize
%     % AB1(n,:) = real(ifft([((AB(n,:))./sqrt(abs(AB(n,:).^2/100+1e10))) conj(fliplr(AB(n,2:end)./sqrt(abs(AB(n,2:end).^2/100+1e10))))]));
%     AB1(:,n) = ((AB(:,n))./sqrt(abs(AB(:,n).^2/100+1e8)));
% end
% % clear AB


% f index for 10 MHz = 110
% f index for 30 MHz = 310
% f index for 74 MHz = 757
% f index for 80 MHz = 820
fc = (0.09765625*110:0.09765625:0.09765625*820)*1e6;
% fc = (0.09765625*600:0.09765625:0.09765625*650)*1e6;
% fc = 0.09765625*390*1e6;
% fc = 60*1e6;

c = physconst('LightSpeed');
% instrument delay
thao = 196 / c;
phasethao = exp(-2*pi*1i*fc'*thao);
AB1 = AB1 .* flip(phasethao);

h = (0:0.083886:0.083886*221723)/3600;
h = h + 5.25;
h = h * pi/12;
% h = (1:-0.083886/3600:-1) * pi/12;

% dec = 0.71094094; % CygnusA Declination
% dec = 1.0262536; % CassA Declination
% dec = [34.200833*pi/180 0.71094094 1.0262536];
% dec = linspace(0.2,1.2,10);
% dec = -10*pi/180;

lat = 34.200833*pi/180; % 34 12 03.000 N LoFASM IV
% lat = 38.433056*pi/180; % 38 25 59.000 N LoFASM III
% lat = 0;
dec = lat + 23*pi/180;
% dec = 10*pi/180;

enu = [0 -152;0 0;0 0];
rotM1 = [0 -sin(lat) cos(lat);1 0 0;0 cos(lat) sin(lat)];
xyz = rotM1 * enu;
lxyz = xyz(:,1) - xyz(:,2);

% visibility2d = zeros(1024,1024,length(dec));
% visibility2dGPU = gpuArray(4000,2000);
uSize = 1024;
vSize = 1024;
visibility = zeros(uSize,vSize);
% visibilitycounter = zeros(uSize,vSize);

for m = 1:length(dec)
    hT = reshape(h,1,1,length(h));
    zeroT = zeros(1,1,length(h));
    decT = dec(m) * ones(1,1,length(h));
    rotM2 = [sin(hT) cos(hT) zeroT;...
            -sin(dec(m))*cos(hT) sin(dec(m))*sin(hT) cos(decT);...
            cos(dec(m))*cos(hT) -cos(dec(m))*sin(hT) sin(decT)];
    uvw = sum(rotM2 .* lxyz',2) .* (fc/c);

    AB2 = reshape(AB1,[1 length(fc)*length(h)]);
    uArray = reshape(uvw(1,:,:),[1 length(fc)*length(h)]);
    vArray = reshape(uvw(2,:,:),[1 length(fc)*length(h)]);
    wArray = reshape(uvw(3,:,:),[1 length(fc)*length(h)]);
    AB2 = AB2 .* exp(-2*pi*1i*wArray);
    
    if (max(uArray) == min(uArray))
        uArray2 = ones(1,length(uArray));
        uaxis = ones(1,uSize) * max(uArray);
    else
        uArray2 = round(uArray/((max(uArray)-min(uArray))/(uSize-1)));
        uArray2 = uArray2 - min(uArray2) + 1;
        uaxis = min(uArray):(max(uArray)-min(uArray))/(uSize-1):max(uArray);
    end
    if (max(vArray) == min(vArray))
        vArray2 = ones(1,length(vArray));
        vaxis = ones(1,vSize) * max(vArray);
    else
        vArray2 = round(vArray/((max(vArray)-min(vArray))/(vSize-1)));
        vArray2 = vArray2 - min(vArray2) + 1;
        vaxis = min(vArray):(max(vArray)-min(vArray))/(vSize-1):max(vArray);
    end
    
    for n = 1:length(uArray2)
        visibility(uArray2(n),vArray2(n)) = visibility(uArray2(n),vArray2(n)) + AB2(n);
        % visibilitycounter(uArray2(n),vArray2(n)) = visibilitycounter(uArray2(n),vArray2(n)) + 1;
    end
    % visibility = visibility ./ visibilitycounter;
    visibility(isnan(visibility)) = 0;
    % visibility2d(:,:,m) = visibility;
end
% visibility(abs(visibility)>10) = 0;
% visibility(abs(visibility)>0.5*1e8) = 0; % when whitening is not used


% % intensity = fftshift(fft2(visibility,512,512));
% intensity = fftshift(fft2(fftshift(visibility)));
% mAxis = linspace(-0.2,0.2,512);
% lAxis = linspace(-pi/4,pi/8,512);
% mesh(lAxis,mAxis,transpose(real(intensity)));


% intensity = zeros(4000,2000);
% for m = 1:length(dec)
%     intensity = intensity + real(fftshift(ifft2(fftshift(visibility2d(:,:,m)))));
% end

% mesh(uaxis,vaxis,abs(transpose(visibility)))
% xlabel('u (wavelength)')
% ylabel('v (wavelength)')
% zlabel('Visibility')
% title('Visibility in u-v Axis')
% visibility(isnan(visibility)) = 0;
% 

%%%% Calculate Full Visibility %%%%
if (abs(min(uArray)) > abs(max(uArray)))
    fulluaxis = min(uArray):(max(uArray)-min(uArray))/(uSize-1):-1*min(uArray);
else
    fulluaxis = -1*max(uArray):(max(uArray)-min(uArray))/(uSize-1):max(uArray);
end
if (abs(min(vArray)) > abs(max(vArray)))
    fullvaxis = min(vArray):(max(vArray)-min(vArray))/(vSize-1):-1*min(vArray);
else
    fullvaxis = -1*max(vArray):(max(vArray)-min(vArray))/(vSize-1):max(vArray);
end

fullVisibility = zeros(length(fulluaxis),length(fullvaxis));

if (abs(min(uArray)) > abs(max(uArray)))
    if (abs(min(vArray)) > abs(max(vArray)))
        fullVisibility(1:length(uaxis),1:length(vaxis)) = visibility;
        fullVisibility(length(fulluaxis)-length(uaxis)+1:end,length(fullvaxis)-length(vaxis)+1:end) = rot90(conj(visibility),2);
    else
        fullVisibility(1:length(uaxis),length(fullvaxis)-length(vaxis)+1:end) = visibility;
        fullVisibility(length(fulluaxis)-length(uaxis)+1:end,1:length(vaxis)) = rot90(conj(visibility),2);
    end    
else
    if (abs(min(vArray)) > abs(max(vArray)))
        fullVisibility(length(fulluaxis)-length(uaxis)+1:end,1:length(vaxis)) = visibility;
        fullVisibility(1:length(uaxis),length(fullvaxis)-length(vaxis)+1:end) = rot90(conj(visibility),2);
    else
        fullVisibility(length(fulluaxis)-length(uaxis)+1:end,length(fullvaxis)-length(vaxis)+1:end) = visibility;
        fullVisibility(1:length(uaxis),1:length(vaxis)) = rot90(conj(visibility),2);
    end
end
% 
% t2 = toc;
% mesh(fulluaxis,fullvaxis,abs(fullVisibility'))
% xlabel('u (wavelength)')
% ylabel('v (wavelength)')
% zlabel('Visibility')
% title('Visibility in u-v Axis')

% %%% Interpolation of Visibility %%%
% [uI,vI] = meshgrid(-20:0.01:20);
% fullVisibilityInterpolated = interp2(fullvaxis,fulluaxis,fullVisibility,uI,vI);
% fullVisibilityInterpolated(isnan(fullVisibilityInterpolated)) = 0;
% % mesh(uI,vI,abs(fullVisibilityInterpolated));
% 
% fulluaxis = uI(1,:);
% fullvaxis = vI(:,1);
% fullVisibility = fullVisibilityInterpolated;
% %%% End of ?nterpolation %%%

nP = 2;
intensityFFTFull = fftshift(fft2(fullVisibility,nP*length(fulluaxis),nP*length(fullvaxis)));
[cU,cV] = size(intensityFFTFull);
ctrU = floor(cU/2);
ctrV = floor(cV/2);
indxU = ctrU-floor(nP*2*max(fulluaxis)):ctrU+floor(nP*2*max(fulluaxis));
indxV = ctrV-floor(nP*2*max(fullvaxis)):ctrV+floor(nP*2*max(fullvaxis));
figure
intensityFFT = intensityFFTFull(indxU,indxV);
[nL,nM] = size(intensityFFT);
mAxis = linspace(-1,1,nM);
lAxis = linspace(-1,1,nL);
mesh(lAxis,mAxis,abs(transpose(intensityFFT)))

% RA = 4.8:(6.2-4.8)/100:6.2;
% Decl = 0.3:(1.2-0.3)/100:1.2;
% intensity = zeros(length(RA),length(Decl));
% for l = 1:length(RA)
%     for m = 1:length(Decl)
%         for u = 1:length(uaxis)
%             for v = 1:length(vaxis)
%                 intensity(l,m) = intensity(l,m) + visibility2(u,v)*exp(1i * 2*pi * (uaxis(u)*l + vaxis(v)*m));
%             end
%         end
%         intensity(l,m) = intensity(l,m)/(length(uaxis)*length(vaxis));
%     end
% end