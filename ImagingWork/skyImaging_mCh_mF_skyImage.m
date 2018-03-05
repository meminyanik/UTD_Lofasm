tStart = 1;

% RFiltered = RFiltered / 1e6;

% AAlT = reshape(RFiltered(1,1,:,:),1024,504);
% ABlT = reshape(RFiltered(2,1,:,:),1024,504);
% BBlT = reshape(RFiltered(2,2,:,:),1024,504);

% AAlT = BBlTerm;
ABlT = ABlTerm / 1024; % * exp(1i*3);
% BBlT = AAlTerm;

% AAlT = AAlT(110:820,tStart:end);
ABlT = ABlT(110:820,tStart:end);
% BBlT = BBlT(110:820,tStart:end);

% AAlT = AAlT(217:417,tStart:end);
% ABlT = ABlT(217:417,tStart:end);
% BBlT = BBlT(217:417,tStart:end);

% ABlT = ABlT(310,tStart:end);

[~, tSize] = size(ABlT);


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
ra = linspace(-20,50,100);

% intensityBF = zeros(length(dec),length(ra));
intensityMF = zeros(length(dec),length(ra));

% fc = (0.09765625*1:0.09765625:0.09765625*1024)'*1e6;
fc = (0.09765625*110:0.09765625:0.09765625*820)'*1e6;
% fc = (0.09765625*217:0.09765625:0.09765625*417)'*1e6;
% fc = (0.09765625*310)'*1e6;
c = physconst('LightSpeed');

% phaseFunctionBF = zeros(length(dec),length(ra),length(fc));
phaseFunctionMF = zeros(length(dec),length(ra),length(fc));

nsF = 67;
thao = nsF/(100e6);
for dD = 1:length(dec)
    for rR = 1:length(ra)

%         tgT = reshape(tgApproximate(dD,rR,:),1,length(h));
        tgT = reshape(tgExact(dD,rR,:,:),6,length(h));
        
        %%% Coherent
%         phaseComp = exp(-2*pi*1i*fc*(tgT-thao));
        phaseComp = 0;
        for nA = 1:6
            phaseComp = phaseComp + exp(-2*pi*1i*fc*(tgT(nA,:)-thao));
        end
        
%         intensityT = AAlT + 2*real(ABlT .* phaseComp) + BBlT;
%         intensityBF(dD,rR) = sum(sum(intensityT));
        
        intensityMF(dD,rR) = sum(sum(ABlT .* phaseComp));
        
%         phaseFunctionBF(dD,rR,:) = sum(intensityT , 2);
        phaseFunctionMF(dD,rR,:) = sum(ABlT .* phaseComp , 2);
        
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
    end
end

mesh(ra,dec,abs(intensityMF),'FaceColor','interp','LineStyle','none')
xlabel('RA (Degrees)')
ylabel('Dec (Degrees)')
colormap('jet'); 

% figure
% mesh(ra,dec,abs(intensityBF))
% xlabel('RA (Degrees)')
% ylabel('Dec (Degrees)')
% colormap('jet'); 
%
figure
CygAPhase = reshape(phaseFunctionMF(26,48,:),1,length(fc));
% plot(fc/1e6,unwrap(angle(CygAPhase)))
plot(fc/1e6,angle(CygAPhase))
xlabel('Frequency (MHz)')
ylabel('Unwrapped Phase Angle (Degrees)')
figure
CasAPhase = reshape(phaseFunctionMF(37,79,:),1,length(fc));
% plot(fc/1e6,unwrap(angle(CasAPhase)))
plot(fc/1e6,angle(CasAPhase))
xlabel('Frequency (MHz)')
ylabel('Unwrapped Phase Angle (Degrees)')