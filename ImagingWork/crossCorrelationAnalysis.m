fSample = 200*1e6; % LoFASM ADCs clock rate is: 200 Msps.
intTime = 0.083886; % integration time (sec);

nSample = fSample * intTime;
t = linspace(0,intTime,nSample);

bandwidth = 0.09765625*1e6; % Hz
fStart = 110;
fStop = 111;

% Design a 20th-order bandpass FIR filter
% The sample rate is fSample Hz.
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',bandwidth*fStart,'CutoffFrequency2',bandwidth*(fStop+1), ...
         'SampleRate',fSample);
% fvtool(bpFilt)

nInt = 120;

analysisSignal = zeros(1024,nInt);


% FM Signal
fCarrier = 90*1e6; % 90 MHz
fDeviation = 200*1e3; % 200 KHz
signal = sin(2*pi*30*t)+2*sin(2*pi*60*t);
fmSignal = fmmod(signal,fCarrier,fSample,fDeviation) + 0.001*randn(1,nSample);

% signalDem = fmdemod(fmSignal,fCarrier,fSample,fDeviation); % Demodulate both channels.
% 
% plot(signal)
% hold on
% plot(signalDem)

% fmSignalFiltered = filter(bpFilt,fmSignal);
% CorrFmSignal = fmSignal .* conj(fmSignal);

for n = 1:nInt
    signalRFI = 0.01 * randn(1,nSample); % + 1i*rand(1,nSample);
    
%     signalRFI = fmSignal + 0.001*randn(1,nSample);
    signalRFI = filter(bpFilt,signalRFI);
    signalRFI = fft(signalRFI,1024) .* conj(fft(signalRFI,1024));

    analysisSignal(n) = sum(signalRFI) / nSample;
end

plot(xcorr(analysisSignal))
