function ABlTerm = longTermAverage(AB,nSAv)
% Obtain Long Term Correlation Matrix
[fSize, tSize] = size(AB);

%%% nSAv: Number of Samples

%%% Use Integration Time
% integrationTime = 10; % 10 Sec
% nSample = round(integrationTime / 0.083886);
% 
% ABpadded = AB;
% lTermtSize = ceil(tSize/nSample);
% ABpadded = [ABpadded zeros(fSize,lTermtSize*nSample-tSize)];
% ABlTerm = zeros(fSize,lTermtSize);
% 
% for k = 1:lTermtSize
%     ABlTerm(:,k) = sum(ABpadded(:,(k-1)*nSample+1:k*nSample),2)/nSample;
% end


%% If there are missing samples
% lTermtSize = floor(tSize/nSAv);
% ABlTerm = zeros(fSize,lTermtSize);
% 
% for k = 1:lTermtSize
%     ABPart = AB(:,(k-1)*nSAv+1:k*nSAv);
%     nonZeroSamples = nSAv - sum(all(ABPart==0,1));
%     if nonZeroSamples == 0
%         ABlTerm(:,k) = 0;
%     else
%         ABlTerm(:,k) = sum(ABPart,2)/nonZeroSamples;
%     end
% end
% 
% % Padding with Mean Value
% ABlTerm(~abs(ABlTerm)) = mean(ABlTerm(:));


%% If there are missing samples
lTermtSize = floor(tSize/nSAv);
ABlTerm = zeros(fSize,lTermtSize);

for k = 1:lTermtSize
    ABPart = AB(:,(k-1)*nSAv+1:k*nSAv);
    ABlTerm(:,k) = sum(ABPart,2)/nSAv;
end