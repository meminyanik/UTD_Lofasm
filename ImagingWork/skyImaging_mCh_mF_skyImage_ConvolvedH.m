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
%-------------------------------------------------------------------------%

% dec = linspace(-pi/2,pi/2,100);
dec = linspace(0,80,50);
% ra = linspace(-pi/2,pi/2,100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For Long Term Correlation
%-------------------------------------------------------------------------%
ra = (raLs*nSAv*0.083886 : nSAv*0.083886 : raRs*nSAv*0.083886)/3600;
ra = ra * 180/12;
%-------------------------------------------------------------------------%

% intensityBF = zeros(length(dec),length(ra));
intensityMF = zeros(length(dec),length(ra));


% fc = (0.09765625*1:0.09765625:0.09765625*1024)'*1e6;
fc = (0.09765625*110:0.09765625:0.09765625*820)'*1e6;
% fc = (0.09765625*217:0.09765625:0.09765625*417)'*1e6;
c = physconst('LightSpeed');

% phaseFunctionBF = zeros(length(dec),length(ra),length(fc));
% phaseFunctionMF = zeros(length(dec),length(ra),length(fc));

[~, padSl] = size(raHl);
[~, padSr] = size(raHr);
ABlT = [zeros(fSize,padSl) ABlT zeros(fSize,padSr)];
ABlT = rot90(ABlT,2);
[fSizeN, tSizeN] = size(ABlT);


%%% Two 1D FFT
% ABlTfft = fft(ABlT,fSizeN,1);
% % ABlTfft = ABlTfft(1:178,:);
% ABlTfft = fft(ABlTfft,tSizeN,2);

%%% 2D FFT
ABlTfft = fft2(ABlT);

nsF = 69.5;
thao = nsF/(100e6);
for dD = 1:length(dec)
    %for rR = 1:length(ra)
       
        % tgT = tgApproximate(dD,rR:rR+length(hT)-1);
        % tgT = reshape(tgExact(dD,:,rR:rR+length(hT)-1),6,length(hT));
        tgT = reshape(tgExact(dD,:,:),6,length(h));

        %%% Coherent
%         phaseComp = exp(-2*pi*1i*fc*(tgT-thao));
        phaseComp = 0;
        for nA = 1:6
            phaseComp = phaseComp + exp(-2*pi*1i*fc*(tgT(nA,:)-thao));
        end
        
        %%% Two 1D FFT
%         phaseCompfft = fft(phaseComp,fSizeN,1);
%         % phaseCompfft = phaseCompfft(1:178,:);
%         phaseCompfft = fft(phaseCompfft,tSizeN,2);
        
        %%% 2D FFT
        phaseCompfft = fft2(phaseComp);
        
%         intensityT = AAlT + 2*real(ABlT .* phaseComp) + BBlT;
%         intensityBF(dD,rR) = sum(sum(intensityT));
        
        intensityMFT = ifft2(ABlTfft.*phaseCompfft);
        intensityMFT = circshift(intensityMFT,padSr,2);
        % intensityMFT = circshift(intensityMFT,1,1);
        intensityMF(dD,:) = intensityMFT(1,1:length(ra));
        % intensityMFT = conv2(ABlT,phaseComp);
        
        % intensityMF(dD,rR) = sum(sum(ABlT .* phaseComp));
        
%         phaseFunctionBF(dD,rR,:) = sum(intensityT , 2);
%         phaseFunctionMF(dD,rR,:) = sum(ABlT .* phaseComp , 2);
        
        
%         %%% Nonceherent
%         phaseComp = exp(-2*pi*1i*fc*tgT);
%         phaseComp = 0;
%         for nA = 1:6
%             phaseComp = phaseComp + exp(-2*pi*1i*fc*tgT(nA,:));
%         end
%         intensityT = AAlT + ABlT .* phaseComp + conj(ABlT .* phaseComp) + BBlT;
%         intensityBF(dD,rR) = sum(abs(sum(intensityT , 2)).^2);
% 
%         intensityMF(dD,rR) = sum(abs(sum(ABlT .* phaseComp , 2)).^2);
    %end
end
intensityMF = fliplr(intensityMF);
mesh(ra,dec,abs(intensityMF))
xlabel('RA (Degrees)')
ylabel('Dec (Degrees)')
colormap('jet'); 

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