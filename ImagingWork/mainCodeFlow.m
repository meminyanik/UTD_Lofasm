load('AB_Padded.mat')
ABwF = whiteningApproach1(AB);
ABwF = ifft(ABwF);
ABlTerm = longTermAverage(ABwF);
ABlTerm = fft(ABlTerm);