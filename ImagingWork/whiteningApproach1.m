function ABwF = whiteningApproach1(AB)
[fSize, tSize] = size(AB);
ABwF = zeros(fSize, tSize);

lambda = 1e8;

%% Original Version
for n = 1:tSize
    ABwF(:,n) = AB(:,n) ./ sqrt(abs(AB(:,n).^2 + lambda));
end

% for n = 1:tSize
%     ABwF(:,n) = AB(:,n) ./ (abs(AB(:,n)).^2 + lambda);
% end