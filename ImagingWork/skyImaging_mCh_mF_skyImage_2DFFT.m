tStart = 1;
load('tgExact.mat')

load('peakPhaseCygA');

% RFiltered = RFiltered / 1e6;

% AAlT = reshape(RFiltered(1,1,:,:),1024,504);
% ABlT = reshape(RFiltered(2,1,:,:),1024,504);
% BBlT = reshape(RFiltered(2,2,:,:),1024,504);

% ABlT = ABSim;

% load('AAlTerm');
% AAlT = AAlTerm;
load('ABlTerm');
ABlT = ABlTerm / sqrt(2048) .* exp(1i*peakPhase)'; % * exp(-1i*2.265);
% load('BBlTerm');
% BBlT = BBlTerm;

% AAlT = AAlT(110:820,tStart:end);
% ABlT = ABlT(110:820,tStart:end);
% BBlT = BBlT(110:820,tStart:end);

% AAlT = AAlT(217:417,tStart:end);
ABlT = ABlT(317:417,tStart:end);
% BBlT = BBlT(217:417,tStart:end);

% ABlT = ABlT(310,tStart:end);

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
h = (0*nSAv*0.083886 : nSAv*0.083886 : (tSize-1)*nSAv*0.083886)/3600;
h = (h + recordStart) * pi/12;
%-------------------------------------------------------------------------%

% dec = linspace(-pi/2,pi/2,100);
dec = linspace(0,80,50);
% ra = linspace(-pi/2,pi/2,100);
ra = linspace(-20,80,100);
dSize = length(dec);
rSize = length(ra);

% intensityBF = zeros(length(dec),length(ra));
intensityMF = zeros(dSize,rSize);

% fc = (0.09765625*1:0.09765625:0.09765625*1024)'*1e6;
% fc = (0.09765625*110:0.09765625:0.09765625*820)'*1e6;
fc = (0.09765625*317:0.09765625:0.09765625*417)'*1e6;
% fc = (0.09765625*310)'*1e6;
c = physconst('LightSpeed');

% phaseFunctionBF = zeros(length(dec),length(ra),length(fc));
% phaseFunctionMF = zeros(length(dec),length(ra),length(fc));

ABlT = reshape(transpose(ABlT),1,tSize,fSize);
fc = reshape(fc,1,1,fSize);

figure
xlabel('RA (Degrees)')
ylabel('Dec (Degrees)')
colormap('jet');

% nsF = 69;
% thao = nsF/(100e6);
thao = 0;

for hH = 1:tSize
    
    %         tgT = reshape(tgApproximate(dD,rR,:),1,length(h));
    tgT = reshape(tgExact(:,:,:,hH),dSize,rSize,6);
    
    %%% Coherent
    %         phaseComp = exp(-2*pi*1i*fc*(tgT-thao));
    phaseComp = 0;
    for nA = 1:6
        tgTA = reshape(tgT(:,:,nA),dSize,rSize);
        phaseComp = phaseComp + exp(-2*pi*1i*fc.*(tgTA-thao));
    end
    ABlTfft = fft2(ABlT(:,hH,:),dSize,rSize);
    intensityMF = intensityMF + sum(ifft2(ABlTfft.*fft2(phaseComp)),3);
    
    %         intensityT = AAlT + 2*real(ABlT .* phaseComp) + BBlT;
    %         intensityBF(dD,rR) = sum(sum(intensityT));
    
    % intensityMF(dD,rR) = sum(sum(ABlT .* phaseComp));
    
    %         phaseFunctionBF(dD,rR,:) = sum(intensityT , 2);
    %         phaseFunctionMF(dD,rR,:) = sum(ABlT .* phaseComp , 2);
    
    %         %%% Nonceherent
    % %         phaseComp = exp(-2*pi*1i*fc*tgT);
    %         phaseComp = 0;
    %         for nA = 1:6
    %             phaseComp = phaseComp + exp(-2*pi*1i*fc*tgT(nA,:));
    %         end
    %         intensityT = AAlT + ABlT .* phaseComp + conj(ABlT .* phaseComp) + BBlT;
    %         intensityBF(dD,rR) = sum(abs(sum(intensityT , 2)).^2);
    %
    %         intensityMF(dD,rR) = sum(abs(sum(ABlT .* phaseComp , 2)).^2);
    
%     mesh(ra,dec,abs(intensityMF))
%     xlabel('RA (Degrees)')
%     ylabel('Dec (Degrees)')
%     colormap('jet');
    
    
    mesh(ra,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
    view(2)
    drawnow;
end
% mesh(ra,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
% view(2)

% figure
% mesh(ra,dec,abs(intensityBF))
% xlabel('RA (Degrees)')
% ylabel('Dec (Degrees)')
% colormap('jet'); 
% 
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