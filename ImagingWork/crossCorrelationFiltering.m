function [ABFiltered, ABMask] = crossCorrelationFiltering(AB)

[fSize,tSize] = size(AB);
ABFiltered = AB;

maxDistance = 1e34;
nSample = 3520;
% nSample = 120;
nBlock = floor(tSize/nSample);

meanAB = abs(mean(mean(AB)));
Signal = xcorr(meanAB*ones(1,nSample),meanAB*ones(1,nSample));

ABMask = zeros(fSize,tSize);

%% Cross Correlation Filtering
for fI = 1:fSize % Each Frequency Index
    for tI = 1:nBlock % Each Time Block Index
        Data = abs(xcorr(AB(fI,(tI-1)*nSample+1:tI*nSample),AB(fI,(tI-1)*nSample+1:tI*nSample)));
        testData = findsignal(Data,Signal,'MaxDistance',maxDistance);
        
        if (isempty(testData))
            ABFiltered(fI,(tI-1)*nSample+1:tI*nSample) = 0;
            ABMask(fI,(tI-1)*nSample+1:tI*nSample) = 1;
        end
    end
end