function PSF = calculatePointSpreadFunction2DFFT(tgExact,decInd,raInd,lengthT,fStart,fStop,nSAv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Frequency Range
%-------------------------------------------------------------------------%
fc = (0.09765625*fStart:0.09765625:0.09765625*fStop)'*1e6;
%-------------------------------------------------------------------------%

nsF = 69.5;
thao = nsF/(100e6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Simulated AB Data - Method 1
%-------------------------------------------------------------------------%
% tgT = reshape(tgExact(decInd,:,raInd:raInd+lengthT-1),6,lengthT);
% phaseComp = 0;
% for nA = 1:6
%     phaseComp = phaseComp + exp(-2*pi*1i*fc*(tgT(nA,:)-thao));
% end 
% ABlT = conj(phaseComp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Simulated AB Data - Method 2
%---------------------------------------------------------------------%
phaseComp = 0;
[~,~,tSizeN] = size(tgExact);
for nA = 1:6
    tgT = reshape(tgExact(decInd,nA,:),1,tSizeN);
    phaseComp = phaseComp + exp(-2*pi*1i*fc.*(tgT-thao));
end
ABlT = conj(phaseComp(:,raInd+1:raInd+lengthT));
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract Related Frequency and Time Range
%-------------------------------------------------------------------------%
ABlT = ABlT(fStart:fStop,:);
[fSize, ~] = size(ABlT);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna Coordinates
%-------------------------------------------------------------------------%
mainLongMap = -116.793272;
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the number of days from J2000.0 for 05:45:11 hrs UT on 2016 July 21th
%-------------------------------------------------------------------------%
startHour = 05 + 45/60 + 11/3600;
fractionDay = startHour / 24;
dayNum = 21;
daysSinceJanuary = 182;
daysSinceJ2000 = 5842.5;
numOfDays = daysSinceJ2000 + daysSinceJanuary + dayNum + fractionDay;

% Find the Local Siderial Time (LST) of Each Antenna
lstMain = 100.46 + 0.985647 * numOfDays + mainLongMap + 15*startHour;
lstMain = mod(lstMain,360);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Total Scanned Hour Duration
%-------------------------------------------------------------------------%
raLd = 0; raRd = 90;

hT = (0*nSAv*0.083886 : nSAv*0.083886 : (lengthT-1)*nSAv*0.083886)/3600;
raLs = round(raLd*12/180*3600/(nSAv*0.083886));
raRs = round(raRd*12/180*3600/(nSAv*0.083886));

raHl = (raLs*nSAv*0.083886 : nSAv*0.083886 : -1*nSAv*0.083886)/3600;
raHr = (lengthT*nSAv*0.083886 : nSAv*0.083886 : (lengthT+raRs)*nSAv*0.083886)/3600;
h = [raHl hT raHr] * pi/12; % in radians
tSizeN = length(h);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Declination (in degrees)
%-------------------------------------------------------------------------%
dec = linspace(0,80,50);
dSize = length(dec);
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define Right Ascention (in hours)
%-------------------------------------------------------------------------%
rA = (raLs*nSAv*0.083886 : nSAv*0.083886 : raRs*nSAv*0.083886)/3600;
rA = rA + lstMain*12/180;
rSize = length(rA);
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare AB for IFFT
%-------------------------------------------------------------------------%
[~, padSl] = size(raHl);
[~, padSr] = size(raHr);
ABlT = [zeros(fSize,padSl) ABlT zeros(fSize,padSr)];
ABlT = fliplr(ABlT);
%-------------------------------------------------------------------------%

ABlT = reshape(transpose(ABlT),1,tSizeN,fSize);
fc = reshape(fc,1,1,fSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D FFT
ABlTfft = fft2(ABlT,dSize,tSizeN);
%-------------------------------------------------------------------------%

nsF = 69.5;
thao = nsF/(100e6);
% thao = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Approximate Matched Filter with 1 Main Antenna is Used
%---------------------------------------------------------------------%
% phaseComp = exp(-2*pi*1i*fc.*(tgApproximate-thao));
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Exact Matched Filter with 6 Main Antenna is Used
%---------------------------------------------------------------------%
phaseComp = 0;
for nA = 1:6
    tgT = reshape(tgExact(:,nA,:),dSize,tSizeN);
    phaseComp = phaseComp + exp(-2*pi*1i*fc.*(tgT-thao));
end
%---------------------------------------------------------------------%

intensityMFT = sum(ifft2(ABlTfft.*fft2(phaseComp)),3);
intensityMF = circshift(intensityMFT,padSl,2);
PSF = intensityMF(:,1:rSize);

% mesh(imag(PSF),'FaceColor','interp','LineStyle','none')
% xlabel('RA Index')
% ylabel('Dec Index')
% view(2)
% colormap('jet');