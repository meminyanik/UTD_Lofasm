function RFiltered = spatialFilteringR(AB,AA,BB,maxEigenvalue,minEigenvalue)
[fSize,tSize] = size(AB);
RFiltered = zeros(2,2,fSize,tSize);

eigThreshold = abs(diff(maxEigenvalue./minEigenvalue));
eigThreshold = [eigThreshold ; zeros(1,tSize)];
for fI = 1:fSize % Each Frequency Index
    for tI = 1:tSize % Each Time Index
        R = [AA(fI,tI) conj(AB(fI,tI));...
             AB(fI,tI) BB(fI,tI)];
        [u, e, ~] = svd(R);
        
%         if (e(1,1)/e(2,2) > 20)
%             RFiltered(:,:,fI,tI) = R - (e(1,1)-e(2,2))*u(:,1)*u(:,1)';
%         else
%             RFiltered(:,:,fI,tI) = R;
%         end

          if (eigThreshold(fI,tI) > 10)
            RFiltered(:,:,fI,tI) = R - (e(1,1)-e(2,2))*u(:,1)*u(:,1)';
          else
            RFiltered(:,:,fI,tI) = R;
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