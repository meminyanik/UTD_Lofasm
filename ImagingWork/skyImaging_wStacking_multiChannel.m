% load('AB_054511_105606.mat')

ABlT = ABlTerm(110:820,:);
[fSize, tSize] = size(ABlT);
lTnSample = 447;

mainLat = 35 + 14/60 + 50.01/3600;
mainLong = 116 + 47/60 + 35.78/3600;
outLat = 35 + 14/60 + 50.93/3600;
outLong = 116 + 47/60 + 29.37/3600;

decCenter = 90*pi/180;

% intensityAbsSum = zeros(length(lq),length(mq));

% f index for 10 MHz = 110
% f index for 30 MHz = 310
% f index for 80 MHz = 820
fc = (0.09765625*110:0.09765625:0.09765625*820)*1e6;
% fc = 60*1e6;
c = physconst('LightSpeed');

recordStart = 5 + 45/60 + 11/3600;
h = (0*lTnSample*0.083886 : lTnSample*0.083886 : (tSize-1)*lTnSample*0.083886)/3600;
% h = (0.083886*0:0.083886:0.083886*225299)/3600;
h = h - ((mainLong/360*24) - recordStart);
h = h * pi/12;
% h = (10*0:10:10*1863)/3600;
% h = (0:0.083886:0.083886*36954)/3600;
% h = h + 5.25;
% h = (1:-0.083886/3600:-1) * pi/12;

% dec = 0.71094094; % CygnusA Declination
% dec = 1.0262536; % CassA Declination
% dec = [34.200833*pi/180 0.71094094 1.0262536];
% dec = linspace(0.2,1.2,10);
% dec = -10*pi/180;

eDistance = 539.84*cosd(10)*0.3048; % East-West baseline in meters
nDistance = 539.84*sind(10)*0.3048; % North-South baseline in meters
enu = [0 eDistance;0 nDistance;0 0];
rotM1 = [0 -sin(mainLat) cos(mainLat);1 0 0;0 cos(mainLat) sin(mainLat)];
xyz = rotM1 * enu;
lxyz = xyz(:,1) - xyz(:,2);

% visibility2d = zeros(1024,1024,length(dec));
% visibility2dGPU = gpuArray(4000,2000);
uSize = 128;
vSize = 128;
wSize = 128;
visibilityWstacking = zeros(uSize,vSize, wSize);
% visibilitycounter = zeros(uSize,vSize);

hT = reshape(h,1,1,length(h));
zeroT = zeros(1,1,length(h));
decT = decCenter * ones(1,1,length(h));
rotM2 = [sin(hT) cos(hT) zeroT;...
    -sin(decCenter)*cos(hT) sin(decCenter)*sin(hT) cos(decT);...
    cos(decCenter)*cos(hT) -cos(decCenter)*sin(hT) sin(decT)];
uvw = sum(rotM2 .* lxyz',2) .* (fc/c);

ABlT = reshape(ABlT,[1 length(fc)*length(h)]);
uArray = reshape(uvw(1,:,:),[1 length(fc)*length(h)]);
vArray = reshape(uvw(2,:,:),[1 length(fc)*length(h)]);
wArray = reshape(uvw(3,:,:),[1 length(fc)*length(h)]);

if (max(uArray) == min(uArray))
    uArray2 = ones(1,length(uArray));
    uAxis = ones(1,uSize) * max(uArray);
else
    uArray2 = round(uArray/((max(uArray)-min(uArray))/(uSize-1)));
    uArray2 = uArray2 - min(uArray2) + 1;
    uAxis = min(uArray):(max(uArray)-min(uArray))/(uSize-1):max(uArray);
end
if (max(vArray) == min(vArray))
    vArray2 = ones(1,length(vArray));
    vAxis = ones(1,vSize) * max(vArray);
else
    vArray2 = round(vArray/((max(vArray)-min(vArray))/(vSize-1)));
    vArray2 = vArray2 - min(vArray2) + 1;
    vAxis = min(vArray):(max(vArray)-min(vArray))/(vSize-1):max(vArray);
end
if (max(wArray) == min(wArray))
    wArray2 = ones(1,length(wArray));
    wAxis = ones(1,wSize) * max(wArray);
else
    wArray2 = round(wArray/((max(wArray)-min(wArray))/(wSize-1)));
    wArray2 = wArray2 - min(wArray2) + 1;
    wAxis = min(wArray):(max(wArray)-min(wArray))/(wSize-1):max(wArray);
end

for n = 1:length(uArray2)
    visibilityWstacking(uArray2(n),vArray2(n),wArray2(n)) = visibilityWstacking(uArray2(n),vArray2(n),wArray2(n)) + ABlT(n);
end
visibilityWstacking(isnan(visibilityWstacking)) = 0;

%%%% Calculate Full u-v-w Axis %%%%
if (abs(min(uArray)) > abs(max(uArray)))
    fullUaxis = min(uArray):(max(uArray)-min(uArray))/(uSize-1):-1*min(uArray);
else
    fullUaxis = -1*max(uArray):(max(uArray)-min(uArray))/(uSize-1):max(uArray);
end
if (abs(min(vArray)) > abs(max(vArray)))
    fullVaxis = min(vArray):(max(vArray)-min(vArray))/(vSize-1):-1*min(vArray);
else
    fullVaxis = -1*max(vArray):(max(vArray)-min(vArray))/(vSize-1):max(vArray);
end
if (abs(min(wArray)) > abs(max(wArray)))
    fullWaxis = min(wArray):(max(wArray)-min(wArray))/(wSize-1):-1*min(wArray);
else
    fullWaxis = -1*max(wArray):(max(wArray)-min(wArray))/(wSize-1):max(wArray);
end

% %%% Full Visibility Calculation Version 1
% fullVisibility = zeros(length(fullUaxis),length(fullVaxis),length(wAxis));
% for vN = 1:length(wAxis)
%     if (abs(min(uArray)) > abs(max(uArray)))
%         if (abs(min(vArray)) > abs(max(vArray)))
%             fullVisibility(1:length(uAxis),1:length(vAxis),vN) = visibilityWstacking(:,:,vN);
%             fullVisibility(length(fullUaxis)-length(uAxis)+1:end,length(fullVaxis)-length(vAxis)+1:end,vN) = rot90(conj(visibilityWstacking(:,:,vN)),2);
%         else
%             fullVisibility(1:length(uAxis),length(fullVaxis)-length(vAxis)+1:end,vN) = visibilityWstacking(:,:,vN);
%             fullVisibility(length(fullUaxis)-length(uAxis)+1:end,1:length(vAxis),vN) = rot90(conj(visibilityWstacking(:,:,vN)),2);
%         end
%     else
%         if (abs(min(vArray)) > abs(max(vArray)))
%             fullVisibility(length(fullUaxis)-length(uAxis)+1:end,1:length(vAxis),vN) = visibilityWstacking(:,:,vN);
%             fullVisibility(1:length(uAxis),length(fullVaxis)-length(vAxis)+1:end,vN) = rot90(conj(visibilityWstacking(:,:,vN)),2);
%         else
%             fullVisibility(length(fullUaxis)-length(uAxis)+1:end,length(fullVaxis)-length(vAxis)+1:end,vN) = visibilityWstacking(:,:,vN);
%             fullVisibility(1:length(uAxis),1:length(vAxis),vN) = rot90(conj(visibilityWstacking(:,:,vN)),2);
%         end
%     end
% end


%%% Full Visibility Calculation Version 2
fullVisibility = zeros(length(fullUaxis),length(fullVaxis),length(fullWaxis));

if (abs(min(wArray)) > abs(max(wArray)))
    if (abs(min(uArray)) > abs(max(uArray)))
        if (abs(min(vArray)) > abs(max(vArray)))
            fullVisibility(1:length(uAxis),1:length(vAxis),1:length(wAxis)) = visibilityWstacking;
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,length(fullVaxis)-length(vAxis)+1:end,length(fullWaxis)-length(wAxis)+1:end) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        else
            fullVisibility(1:length(uAxis),length(fullVaxis)-length(vAxis)+1:end,1:length(wAxis)) = visibilityWstacking;
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,1:length(vAxis),length(fullWaxis)-length(wAxis)+1:end) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        end
    else
        if (abs(min(vArray)) > abs(max(vArray)))
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,1:length(vAxis),1:length(wAxis)) = visibilityWstacking;
            fullVisibility(1:length(uAxis),length(fullVaxis)-length(vAxis)+1:end,length(fullWaxis)-length(wAxis)+1:end) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        else
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,length(fullVaxis)-length(vAxis)+1:end,1:length(wAxis)) = visibilityWstacking;
            fullVisibility(1:length(uAxis),1:length(vAxis),length(fullWaxis)-length(wAxis)+1:end) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        end
    end
else
    if (abs(min(uArray)) > abs(max(uArray)))
        if (abs(min(vArray)) > abs(max(vArray)))
            fullVisibility(1:length(uAxis),1:length(vAxis),length(fullWaxis)-length(wAxis)+1:end) = visibilityWstacking;
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,length(fullVaxis)-length(vAxis)+1:end,1:length(wAxis)) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        else
            fullVisibility(1:length(uAxis),length(fullVaxis)-length(vAxis)+1:end,length(fullWaxis)-length(wAxis)+1:end) = visibilityWstacking;
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,1:length(vAxis),1:length(wAxis)) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        end
    else
        if (abs(min(vArray)) > abs(max(vArray)))
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,1:length(vAxis),length(fullWaxis)-length(wAxis)+1:end) = visibilityWstacking;
            fullVisibility(1:length(uAxis),length(fullVaxis)-length(vAxis)+1:end,1:length(wAxis)) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        else
            fullVisibility(length(fullUaxis)-length(uAxis)+1:end,length(fullVaxis)-length(vAxis)+1:end,length(fullWaxis)-length(wAxis)+1:end) = visibilityWstacking;
            fullVisibility(1:length(uAxis),1:length(vAxis),1:length(wAxis)) = flip(flip(flip(conj(visibilityWstacking),1),2),3);
        end
    end
end


[lq,mq] = meshgrid(-1:0.005:1);
intensityPerW = zeros(length(lq),length(mq),length(fullWaxis));

for nV = 1:length(fullWaxis)
    nP = 2;
    intensityFFTFull = fftshift(ifft2(fullVisibility(:,:,nV),nP*length(fullUaxis),nP*length(fullVaxis)));
    [cU,cV] = size(intensityFFTFull);
    ctrU = floor(cU/2);
    ctrV = floor(cV/2);
    indxU = ctrU-floor(nP*2*max(fullUaxis)):ctrU+floor(nP*2*max(fullUaxis));
    indxV = ctrV-floor(nP*2*max(fullVaxis)):ctrV+floor(nP*2*max(fullVaxis));
    intensityFFT = intensityFFTFull(indxU,indxV);
    [nL,nM] = size(intensityFFT);
    mAxis = linspace(-1,1,nM);
    lAxis = linspace(-1,1,nL);
    [l,m] = meshgrid(lAxis,mAxis);
    % intensityAbs = abs(intensityFFT');
    intensityReal = real(intensityFFT');
    intensityImag = imag(intensityFFT');
    % intensityInterpolatedAbs = interp2(l,m,intensityAbs,lq,mq);
    intensityInterpolatedReal = interp2(l,m,intensityReal,lq,mq);
    intensityInterpolatedImag = interp2(l,m,intensityImag,lq,mq);
    intensityInterpolatedComplex = intensityInterpolatedReal + 1i*intensityInterpolatedImag;
    
    intensityPerW(:,:,nV) = intensityInterpolatedComplex;
    %     intensityAbsSum = intensityAbsSum + intensityInterpolatedAbs;
    %     % mesh(asin(lq),asin(mq),intensity)
    %     mesh(lq,mq,fliplr(intensityAbsSum))
    %     figure
    %     % mesh(asin(lq),asin(mq),sum(abs(intensity3d),3))
    %     % mesh(lq,mq,sum(abs(intensityPerFrequency),3))
    %     mesh(lq,mq,fliplr(abs(sum(intensityPerW,3))))
    %     t2 = toc;
end

l = -1:0.005:1;
m = -1:0.005:1;
intensityPerWT = intensityPerW;
dirCosN = sqrt(1 - lq.^2 - mq.^2);
for lN = 1:length(l)
    for mN = 1:length(m)
        if (l(lN)^2 + m(mN)^2) > 1
            dirCosN(lN,mN) = 0;
            intensityPerWT(lN,mN,:) = 0;
        end
    end
end

wAxisReshaped = reshape(fullWaxis,1,1,length(fullWaxis));
phaseComp = exp(2*pi*1i*wAxisReshaped .* (dirCosN-1));

intensityPerWT = intensityPerWT .* phaseComp;
intensityT = sum(intensityPerWT,3);
intensity = intensityT .* dirCosN / (max(fullWaxis) - min(fullWaxis));

mesh(l,m,abs(intensity))

% intensityCoherent = zeros(length(lq),length(mq));
% % intensity3d2 = zeros(length(lq),length(mq),fSize);
% for nF = 110:820
%     fc = 0.09765625*nF*1e6;
%     nsF = 75;
%     thao = 1/(0.09765625*nsF*1e6);
%     % thao = 1.6 * 30 / c;
%     intensityCoherent = intensityCoherent + intensityPerFrequency(:,:,nF) * exp(2*pi*1i*fc*thao);
%     % intensity3d2(:,:,nF) = intensity3d(:,:,nF) * exp(2*pi*1i*fc*thao);
%     % intensity2 = intensity2 + intensity3d(:,:,nF) * conj(intensityPhase(nF));
%     % intensity2 = intensity2 + abs(intensity3d(:,:,nF));
%     % mesh(asin(lq),asin(mq),abs(intensity3d(:,:,nF)))
% end
% figure
% % mesh(asin(lq),asin(mq),abs(intensity2))
% mesh(lq,mq,fliplr(abs(intensityCoherent)))
% % figure
% % mesh(lq,mq,abs(sum(intensity3d2,3)))