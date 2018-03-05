tic
% % % load('AB_054511_084543.mat')
% load('AB_054511_105606.mat')
% t1 = toc;

% ABf = smooth(AB,'rlowess');
% ABf = reshape(ABf,fSize, tSize);
ABperFrequency = ABwF(110,:);
[fSize, tSize] = size(ABperFrequency);

[lq,mq] = meshgrid(-1:0.005:1);
intensityPerFrequency = zeros(length(lq),length(mq),fSize);
intensityAbsSum = zeros(length(lq),length(mq));

ABt = ABwF;

for nF = 110:820
    AB1 = ABt(nF,:);
%     varAB1 = var(AB1);
%     AB1(AB1>2*sqrt(varAB1))=0;
%     AB1 = smooth(AB1);
    
    % f index for 10 MHz = 110
    % f index for 30 MHz = 310
    % f index for 80 MHz = 820
    % fc = (0.09765625*110:0.09765625:0.09765625*820)*1e6;
    % fc = (0.09765625*600:0.09765625:0.09765625*650)*1e6;
    fc = 0.09765625*nF*1e6;
    % fc = 60*1e6;
    % h = (0:0.083886:0.083886*221723)/3600;
    h = (0.083886*0:0.083886:0.083886*221723)/3600;
    % h = (10*0:10:10*1863)/3600;
    % h = (0:0.083886:0.083886*36954)/3600;
    % h = h + 5.25;
    h = h * pi/12;
    % h = (1:-0.083886/3600:-1) * pi/12;
    c = physconst('LightSpeed');
    
    % dec = 0.71094094; % CygnusA Declination
    % dec = 1.0262536; % CassA Declination
    % dec = [34.200833*pi/180 0.71094094 1.0262536];
    % dec = linspace(0.2,1.2,10);
    % dec = -10*pi/180;
    
    lat = 34.200833*pi/180; % LoFASM IV   34? 12’ 03.000” N       118? 10’ 18.000” W    -8.35? 
    % lat = 38.433056*pi/180; % LoFASM III   38? 25’ 59.000” N       79? 50’ 23.000” W     -19.7?
    % lat = 0;
    dec = 90*pi/180;
    % dec = 10*pi/180;
    
    enu = [0 -152;0 0;0 0];
    rotM1 = [0 -sin(lat) cos(lat);1 0 0;0 cos(lat) sin(lat)];
    xyz = rotM1 * enu;
    lxyz = xyz(:,1) - xyz(:,2);
    
    % visibility2d = zeros(1024,1024,length(dec));
    % visibility2dGPU = gpuArray(4000,2000);
    uSize = 128;
    vSize = 128;
    % wSize = 128;
    visibility = zeros(uSize,vSize);
    % visibilityWstacking = zeros(uSize,vSize, wSize);
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

%           uArray = cos(h) * 152 * (fc/c);
%           vArray = sin(dec(m)) * sin(h) * 152 * (fc/c);
%           wArray = -cos(dec(m)) * sin(h) * 152 * (fc/c);

%         AB2 = AB2 .* exp(2*pi*1i*wArray);
%         thao = 7.5 * 30 / c;
%         AB2 = AB2 * exp(2*pi*1i*fc*thao);
        
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
%         if (max(wArray) == min(wArray))
%             wArray2 = ones(1,length(wArray));
%             waxis = ones(1,wSize) * max(wArray);
%         else
%             wArray2 = round(wArray/((max(wArray)-min(wArray))/(wSize-1)));
%             wArray2 = wArray2 - min(wArray2) + 1;
%             waxis = min(wArray):(max(wArray)-min(wArray))/(wSize-1):max(wArray);
%         end
        
        for n = 1:length(uArray2)
            visibility(uArray2(n),vArray2(n)) = visibility(uArray2(n),vArray2(n)) + AB2(n);
            % visibilityWstacking(uArray2(n),vArray2(n),wArray2(n)) = visibilityWstacking(uArray2(n),vArray2(n),wArray2(n)) + AB2(n);
            % visibilitycounter(uArray2(n),vArray2(n)) = visibilitycounter(uArray2(n),vArray2(n)) + 1;
        end
        % visibility = visibility ./ visibilitycounter;
        visibility(isnan(visibility)) = 0;
        % visibilityWstacking(isnan(visibilityWstacking)) = 0;
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
    % %%% End of Interpolation %%%
    
    nP = 2;
    intensityFFTFull = fftshift(ifft2(fullVisibility,nP*length(fulluaxis),nP*length(fullvaxis)));
%     intensityFFTFull = fftshift(ifft2(visibility,nP*length(uaxis),nP*length(vaxis)));
    [cU,cV] = size(intensityFFTFull);
    ctrU = floor(cU/2);
    ctrV = floor(cV/2);
    indxU = ctrU-floor(nP*2*max(fulluaxis)):ctrU+floor(nP*2*max(fulluaxis));
    indxV = ctrV-floor(nP*2*max(fullvaxis)):ctrV+floor(nP*2*max(fullvaxis));
%     indxU = ctrU-floor(nP*2*min(uaxis)):ctrU+floor(nP*2*max(uaxis));
%     indxV = ctrV-floor(nP*2*min(vaxis)):ctrV+floor(nP*2*max(vaxis));
    % figure
    intensityFFT = intensityFFTFull(indxU,indxV);
    [nL,nM] = size(intensityFFT);
    mAxis = linspace(-1,1,nM);
    lAxis = linspace(-1,1,nL);
    [l,m] = meshgrid(lAxis,mAxis);
    intensityAbs = abs(intensityFFT');
    intensityReal = real(intensityFFT');
    intensityImag = imag(intensityFFT');
    intensityInterpolatedAbs = interp2(l,m,intensityAbs,lq,mq);
    intensityInterpolatedReal = interp2(l,m,intensityReal,lq,mq);
    intensityInterpolatedImag = interp2(l,m,intensityImag,lq,mq);
    intensityInterpolatedComplex = intensityInterpolatedReal + 1i*intensityInterpolatedImag;
    % figure
    % mesh(lq,mq,intensityInterpolated)
    
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
    
    intensityPerFrequency(:,:,nF) = intensityInterpolatedComplex;
    intensityAbsSum = intensityAbsSum + intensityInterpolatedAbs;
end
% mesh(asin(lq),asin(mq),intensity)
mesh(lq,mq,fliplr(intensityAbsSum))
figure
% mesh(asin(lq),asin(mq),sum(abs(intensity3d),3))
% mesh(lq,mq,sum(abs(intensityPerFrequency),3))
mesh(lq,mq,fliplr(abs(sum(intensityPerFrequency,3))))
t2 = toc;


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