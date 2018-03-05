deltaF = 0.09765625*1e6; % 100e6/1024
indF = 0:1023;

tI = 60/100e6;
x = exp(1i*2*pi*indF*deltaF*tI);

nFFT = 1024;
sfft = abs(fftshift(fft(x,nFFT)));
plot(sfft)
[~,ind] = max(sfft);
tIEst = (ind-1-nFFT/2)/(nFFT*deltaF);

% n = length(indF);
% sfftM = zeros(1,nFFT);
% for k = 1:nFFT
%     for j = 1:n
%         sfftM(k) = sfftM(k) + x(j)*exp(-1i*2*pi/n*(j-1)*(k-1));
%     end
% end
% figure
% plot(abs(sfftM))

