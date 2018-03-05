function [AB,ABlTerm,ABwF2] = whiteningApproach2_2
AB = [];
ABwF2 = [];
ABlTerm = [];
D = dir('..\lofasm4_outrigger_AB');
for k = 1:62
    ABtemp = transpose(csvread(['..\lofasm4_outrigger_AB\',D(k+2).name]));
    [~, tSize] = size(ABtemp);
    ABtempNorm = sqrt(sum(abs(ABtemp).^2,2)/tSize);
    ABlTerm = [ABlTerm sum(ABtemp,2)];
    AB = [AB ABtemp];
    ABwF2 = [ABwF2 ABtemp./ABtempNorm];
end