function [maxEigenvalue, minEigenvalue] = spatialFilteringEigenvalueAnalysis(AB,AA,BB)
[fSize,tSize] = size(AB);

maxEigenvalue = zeros(fSize,tSize);
minEigenvalue = zeros(fSize,tSize);

for fI = 1:fSize % Each Frequency Index
    for tI = 1:tSize % Each Time Index
        R = [AA(fI,tI) conj(AB(fI,tI));...
             AB(fI,tI) BB(fI,tI)];
        e = svd(R);
        maxEigenvalue(fI,tI) = e(1);
        minEigenvalue(fI,tI) = e(2);
    end
end

% %% Directly From R
% for fI = 1:fSize % Each Frequency Index
%     for tI = 1:tSize % Each Time Index
%         e = svd(R(:,:,fI,tI));
%         maxEigenvalue(fI,tI) = e(1);
%         minEigenvalue(fI,tI) = e(2);
%     end
% end

%%% fIndex = 0:100/1023:100;
% fIndex = linspace(0,100,1024);
%%% tIndex = 0:0.083886:0.083886*3575;
% tIndex = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
% mesh(fIndex,tIndex,mag2db(maxEigenvalueFiltered+0.0001))
% xlabel('Frequency (MHz)')
% ylabel('Time (second)')
% zlabel('Power (dB)')
% title('2D Spectrogram of the Maximum Eigenvalue After Filtering')