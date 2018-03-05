function [AAFiltered, BBFiltered, ABFiltered, ABMask] = spatialFiltering(AA,BB,AB)
[fSize,tSize] = size(AB);
AAFiltered = zeros(fSize,tSize);
BBFiltered = zeros(fSize,tSize);
ABFiltered = zeros(fSize,tSize);

ABMask = zeros(fSize,tSize);

%% Spatial Filtering
for fI = 1:fSize % Each Frequency Index
    for tI = 1:tSize % Each Time Index
        R = [AA(fI,tI) conj(AB(fI,tI));...
             AB(fI,tI) BB(fI,tI)];
        [u, e, ~] = svd(R);
        
        if (e(1,1)/e(2,2) > 2)
            %RFiltered = R - (e(1,1)-e(2,2))*u(:,1)*u(:,1)';
            RFiltered = R - e(1,1)*u(:,1)*u(:,1)';
            AAFiltered(fI,tI) = RFiltered(1,1);
            BBFiltered(fI,tI) = RFiltered(2,2);
            ABFiltered(fI,tI) = RFiltered(2,1);
            ABMask(fI,tI) = 1;
        else
            AAFiltered(fI,tI) = R(1,1);
            BBFiltered(fI,tI) = R(2,2);
            ABFiltered(fI,tI) = R(2,1);
        end
    end
end

%%% fIndex = 0:100/1023:100;
% fIndex = linspace(0,100,1024);
%%% tIndex = 0:0.083886:0.083886*3575;
% tIndex = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
% mesh(fIndex,tIndex,mag2db(maxEigenvalueFiltered+0.0001))
% xlabel('Frequency (MHz)')
% ylabel('Time (second)')
% zlabel('Power (dB)')
% title('2D Spectrogram of the Maximum Eigenvalue After Filtering')