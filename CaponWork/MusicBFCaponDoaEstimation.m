nElements = 10; % Number of Elements
d = 0.5; % Distance between element in wavelength

Ang = [30,40,55]; % DoA of signals between 0-90 and in degrees
nSignal = length(Ang); % number of simulated signals

SNR = 0; % SNR in dB

steeringVectorMatrix = exp(-1i*2*pi*d*(0:nElements-1)'*sind(Ang));

fc = [600;800;1000]; % carrier frequency of simulated signal in Hz
fs = 8000; % sample frequency in Hz
t = linspace(0,1,fs);
signal = exp(1i*2*pi*fc*t);
nSample = length(signal);

% noise = 10^(-SNR/20) * 1/sqrt(2) * (randn(nElements,nSample)+1i*randn(nElements,nSample));
% receivedSignal = steeringVectorMatrix * signal + noise;
receivedSignal = awgn(steeringVectorMatrix*signal,SNR);

R = receivedSignal*receivedSignal'/nSample;
% [Q,D] = eig(R);
% [D,I] = sort(diag(D),1,'descend');
% Q = Q (:,I);
% Qs = Q (:,1:nSignal);
% Qn = Q(:,nSignal+1:nElements);
[Q,D] = eig(R,'vector');
Qn = Q(:,1:nElements-nSignal);

scanAngles = (-90:0.1:90);
a = exp(-1i*2*pi*d*(0:nElements-1)'*sind(scanAngles));
a = a./vecnorm(a);

pMUSIC = zeros(1,length(scanAngles));
pBEAM = zeros(1,length(scanAngles));
pCAPON = zeros(1,length(scanAngles));
for k=1:length(scanAngles)
%Compute MUSIC, Beamforming and MVDR “spectrum”
    pMUSIC(k) = 1 / (a(:,k)'*Qn*Qn'*a(:,k));
    pBEAM(k) = (a(:,k)'*R*a(:,k));
    pCAPON(k) = 1 / (a(:,k)'*inv(R)*a(:,k));
    
    % pMUSIC(k) = (a1(:,k)'*a1(:,k)) / (a1(:,k)'*Qn*Qn'*a1(:,k));
    % pBEAM(k) = (a1(:,k)'*R*a1(:,k)) / (a1(:,k)'*a1(:,k));
end
plot(scanAngles,mag2db(abs(pMUSIC)))
figure
plot(scanAngles,mag2db(abs(pBEAM)))
figure
plot(scanAngles,mag2db(abs(pCAPON)));

%%% Analysis Section
R_inv = inv(R);

% real(det(R))
% prod(D)

% [~,D_inv] = eig(R_inv,'vector');
% 1./D
% D_inv

% min(D)
% min(pBEAM)
% max(D)
% max(pBEAM)

% pCAPON_inv = 1./pCAPON;
% min(D_inv)
% min(pCAPON_inv)
% max(D_inv)
% max(pCAPON_inv)

R_inv_diag = eye(nElements,nElements).*R_inv;
R_inv_offDiag = ~eye(nElements,nElements).*R_inv;

R_diag = eye(nElements,nElements).*R;
R_offDiag = ~eye(nElements,nElements).*R;

pCAPON_approx = zeros(1,length(scanAngles));
pCAPON_test = zeros(1,length(scanAngles));
for k=1:length(scanAngles)
%Compute Approximate MVDR “spectrum”
    % pCAPON_approx(k) = 1 / (a(:,k)'*(R_inv_diag+R_inv_offDiag)*a(:,k));
    
    % pCAPON_approx(k) = 1 / ...
    %    ((a(:,k)'*R_inv_diag*a(:,k)) * (1 + (a(:,k)'*R_inv_offDiag*a(:,k))/(a(:,k)'*R_inv_diag*a(:,k))));
    
    pCAPON_test(k) = (a(:,k)'*-R_inv_offDiag*a(:,k))/(a(:,k)'*R_inv_diag*a(:,k));
    % pCAPON_approx(k) = 1 / ((a(:,k)'*R_inv_diag*a(:,k)) * (1 - pCAPON_test(k)));
    pCAPON_approx(k) = 1 / (a(:,k)'*R_inv_diag*a(:,k)) * (1 + pCAPON_test(k));
end
figure
plot(scanAngles,mag2db(abs(pCAPON_approx)));

[pks,locs] = findpeaks(abs(pCAPON_approx),'SortStr','descend');
estimatedAngles = (scanAngles(locs(1:nSignal)));







% receivedSignalFFT = [receivedSignal ; zeros((length(angles)-nElements),nSample)];
% receivedSignalPERIODOGRAM = [zeros((length(angles)-nElements),nSample) ; receivedSignal];
% pPERIODOGRAM = periodogram(receivedSignalPERIODOGRAM,[],length(angles));
% pFFT = fft(R,length(angles));
% pFFT = sum(pFFT');
% figure
% plot(angles,abs(pFFT));
% figure
% plot(angles,abs(pPERIODOGRAM));