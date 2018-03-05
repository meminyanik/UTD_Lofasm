function [ABFiltered, ABMask, ABCoeff] = crossCorrelationFilteringV2(AB)

[fSize,tSize] = size(AB);
ABFiltered = AB;

nSample = 3520;
% nSample = 120;
nBlock = floor(tSize/nSample);

ABMask = zeros(fSize,tSize);
ABCoeff = zeros(fSize,tSize);

%% Cross Correlation Filtering
for fI = 1:fSize % Each Frequency Index
    for tI = 1:nBlock % Each Time Block Index
        Data = abs(xcorr(AB(fI,(tI-1)*nSample+1:tI*nSample),AB(fI,(tI-1)*nSample+1:tI*nSample)));
        
        ABCoeff(fI,(tI-1)*nSample+1:tI*nSample) = Data(nSample)/Data(nSample+50);
        
        if ((Data(nSample)/Data(nSample+10)) > 4)
            ABFiltered(fI,(tI-1)*nSample+1:tI*nSample) = 0;
            ABMask(fI,(tI-1)*nSample+1:tI*nSample) = 1;
        end
    end
end