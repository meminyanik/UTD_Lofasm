tStart = 1;
load('tgExactConvolvedH.mat')

% RFiltered = RFiltered / 1e6;

% AAlT = reshape(RFiltered(1,1,:,:),1024,504);
% ABlT = reshape(RFiltered(2,1,:,:),1024,504);
% BBlT = reshape(RFiltered(2,2,:,:),1024,504);

ABlT = ABSim;

% load('AAlTerm');
% AAlT = AAlTerm;
% load('ABlTerm');
% ABlTerm = whiteningApproach1(ABlTerm);
% ABlT = ABlTerm / 1024;
% load('BBlTerm');
% BBlT = BBlTerm;

% AAlT = BBlTerm;
% ABlT = ABlTerm / 1024;% * exp(-1i*2.265);
% BBlT = AAlTerm;

% AAlT = AAlT(110:820,tStart:end);
% ABlT = ABlT(110:820,tStart:end);
% ABlT = ABlT.* CygAPhase';
% BBlT = BBlT(110:820,tStart:end);

% ABlT = ABlT(310,tStart:end);

% AAlT = AAlT(217:417,tStart:end);
% ABlT = ABlT(217:417,tStart:end);
% BBlT = BBlT(217:417,tStart:end);

[fSize, tSize] = size(ABlT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Should Be Set - Number of Samples Used During Averaging
%-------------------------------------------------------------------------%
nSAv = 447;
% nSAv = 50;
%-------------------------------------------------------------------------%

recordStart = (5 + 45/60 + 11/3600);
% recordStart = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For Long Term Correlation
%-------------------------------------------------------------------------%
hT = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
raLd = -20; raRd = 50;
raLs = round(raLd*12/180*3600/(nSAv*0.083886));
raRs = round(raRd*12/180*3600/(nSAv*0.083886));

raHl = (raLs*nSAv*0.083886 : nSAv*0.083886 : -1*nSAv*0.083886)/3600;
raHr = (tSize*nSAv*0.083886 : nSAv*0.083886 : (tSize+raRs)*nSAv*0.083886)/3600;
h = [raHl hT raHr];
h = (h + recordStart) * pi/12;
tSizeN = length(h);
%-------------------------------------------------------------------------%

% dec = linspace(-pi/2,pi/2,100);
dec = linspace(0,80,50);
% ra = linspace(-pi/2,pi/2,100);
dSize = length(dec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For Long Term Correlation
%-------------------------------------------------------------------------%
ra = (raLs*nSAv*0.083886 : nSAv*0.083886 : raRs*nSAv*0.083886)/3600;
ra = ra * 180/12;
rSize = length(ra);
%-------------------------------------------------------------------------%

fc = (0.09765625*1:0.09765625:0.09765625*1024)'*1e6;
% fc = (0.09765625*110:0.09765625:0.09765625*820)'*1e6;
% fc = (0.09765625*217:0.09765625:0.09765625*417)'*1e6;
% fc = (0.09765625*311)'*1e6;
c = physconst('LightSpeed');

% phaseFunctionBF = zeros(length(dec),length(ra),length(fc));
% phaseFunctionMF = zeros(length(dec),length(ra),length(fc));

[~, padSl] = size(raHl);
[~, padSr] = size(raHr);
ABlT = [zeros(fSize,padSl) ABlT zeros(fSize,padSr)];
% ABlT = rot90(ABlT,2);
ABlT = fliplr(ABlT);
% [~, tSizeN] = size(ABlT);

ABlT = reshape(transpose(ABlT),1,tSizeN,fSize);
fc = reshape(fc,1,1,fSize);

%%% Two 1D FFT
% ABlTfft = fft(ABlT,fSizeN,1);
% ABlTfft = ABlTfft(1:178,:);
% ABlTfft = fft(ABlTfft,tSizeN,2);

%%% 2D FFT
ABlTfft = fft2(ABlT,dSize,tSizeN);

nsF = 69.5;
thao = nsF/(100e6);
% thao = 0;
    
%%% 1 Approximate Antenna
% phaseComp = exp(-2*pi*1i*fc.*(tgApproximate-thao));

%%% 6 Exact Antennas
phaseComp = 0;
for nA = 1:6
    tgT = reshape(tgExact(:,nA,:),dSize,tSizeN);
    phaseComp = phaseComp + exp(-2*pi*1i*fc.*(tgT-thao));
end

intensityMFT = sum(ifft2(ABlTfft.*fft2(phaseComp)),3);

intensityMF = circshift(intensityMFT,padSr,2);
% intensityMF = circshift(intensityMF,1,1);
intensityMF = intensityMF(:,1:rSize);

% mesh(ra,dec,real(fliplr(intensityMF)),'FaceColor','interp','LineStyle','none')
mesh(ra,dec,real(intensityMF),'FaceColor','interp','LineStyle','none')
xlabel('RA (Degrees)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');

% figure
% mesh(abs(intensityMFT),'FaceColor','interp','LineStyle','none')
% view(2)
% colormap('jet');

% mesh(h,dec,abs(intensityMFT))

% %%% 1D FFT in Time
% % ABlTfft = fft(ABlT,tSizeN,2);
% 
% for fF = 1:length(fc)
%     %for rR = 1:length(ra)
%        
%         % tgT = tgApproximate(dD,rR:rR+length(hT)-1);
%         % tgT = reshape(tgExact(dD,:,rR:rR+length(hT)-1),6,length(hT));
%         % tgT = reshape(tgExact(dD,:,:),6,length(h));
% 
%         %%% Coherent
% %         phaseComp = exp(-2*pi*1i*fc*(tgT-thao));
%         
%         
%         %%% Two 1D FFT
% %         phaseCompfft = fft(phaseComp,fSizeN,1);
% %         % phaseCompfft = phaseCompfft(1:178,:);
% %         phaseCompfft = fft(phaseCompfft,tSizeN,2);
%         
%         %%% 2D FFT
%         phaseCompfft = fft2(phaseComp);
% 
%         %%% 1D FFT in Time
%         % phaseCompfft = fft(phaseComp,tSizeN,2);
%         
% %         intensityT = AAlT + 2*real(ABlT .* phaseComp) + BBlT;
% %         intensityBF(dD,rR) = sum(sum(intensityT));
%         
%         intensityMFT = ABlTfft(fF,:).*conj(phaseCompfft);
% %         intensityMFT = fft(intensityMFT,length(dec),1);
% %         
% %         intensityMF = intensityMF + abs(ifft(intensityMFT,tSizeN,2));
%         intensityMF = intensityMF + ifft2(intensityMFT);
%         
%         
%         % intensityMFT = conv2(ABlT,phaseComp);
%         
%         % intensityMF(dD,rR) = sum(sum(ABlT .* phaseComp));
%         
% %         phaseFunctionBF(dD,rR,:) = sum(intensityT , 2);
%         % phaseFunctionMF(dD,rR,:) = sum(ABlT .* phaseComp , 2);
%         
%         
% %         %%% Nonceherent
% %         phaseComp = exp(-2*pi*1i*fc*tgT);
% %         phaseComp = 0;
% %         for nA = 1:6
% %             phaseComp = phaseComp + exp(-2*pi*1i*fc*tgT(nA,:));
% %         end
% %         intensityT = AAlT + ABlT .* phaseComp + conj(ABlT .* phaseComp) + BBlT;
% %         intensityBF(dD,rR) = sum(abs(sum(intensityT , 2)).^2);
% % 
% %         intensityMF(dD,rR) = sum(abs(sum(ABlT .* phaseComp , 2)).^2);
%     %end
% end
% % intensityMF = fliplr(intensityMF);
% % mesh(ra,dec,abs(intensityMF))
% % mesh(abs(ifft2(intensityMF)))
% mesh(abs(intensityMF))
% xlabel('RA (Degrees)')
% ylabel('Dec (Degrees)')
% colormap('jet'); 

% figure
% mesh(ra,dec,abs(intensityBF))
% xlabel('RA (Degrees)')
% ylabel('Dec (Degrees)')
% colormap('jet'); 

% CygAPhase = reshape(phaseFunctionMF(26,48,:),1,length(fc));
% figure
% plot(fc/1e6,unwrap(angle(CygAPhase)))
% xlabel('Frequency (MHz)')
% ylabel('Unwrapped Phase Angle (Degrees)')
% figure
% CasAPhase = reshape(phaseFunctionMF(38,77,:),1,length(fc));
% plot(fc/1e6,unwrap(angle(CasAPhase)))
% xlabel('Frequency (MHz)')
% ylabel('Unwrapped Phase Angle (Degrees)')