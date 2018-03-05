function ABwF = whiteningApproach2(AB)

[fSize, tSize] = size(AB);
integrationTime = 10; % 10 Sec
nSample = round(integrationTime / 0.083886);

ABpadded = AB;
lTermtSize = ceil(tSize/nSample);
ABpadded = [ABpadded zeros(fSize,lTermtSize*nSample-tSize)];
[fSizePadded, tSizePadded] = size(ABpadded);
ABwF = zeros(fSizePadded,tSizePadded);

for k = 1:lTermtSize
    ABtemp = ABpadded(:,(k-1)*nSample+1:k*nSample);
    [~, tSizeTemp] = size(ABtemp);
    ABtempNorm = sqrt(sum(abs(ABtemp).^2,2)/tSizeTemp);
    ABwF(:,(k-1)*nSample+1:k*nSample) = ABtemp./ABtempNorm;
end

ABwF = ABwF(:,1:tSize);