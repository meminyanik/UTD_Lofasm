% find_max: The Hogbom algoritm begins by finding the highest intensity point
% (in terms of absolute value) on the dirty image.
% Our implementation find the maximum, and then the array index of the maximum
% and record the maximum value at this position on the image model.

% sub_beam: In each iteration, the dirty beam is scaled, positioned at the
% maximum point found in step 1, and subtracted from the dirty image.

% add_noise

intensityDirty = intensityMF;
% intensityDirty = intensityBF;
intensityClean = intensityDirty;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When 2DFFT Imaging is used
%---------------------------------------------------------------------%
indCygDec = 26; indCygRa = 163;
indCasDec = 38; indCasRa = 478;
nSources = 2;
indSources = [indCygDec, indCygRa ; indCasDec, indCasRa];
%---------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% When Standard Imaging is used
%---------------------------------------------------------------------%
% indCygDec = 26; indCygRa = 29;
% indCasDec = 37; indCasRa = 83;
% nSources = 2;
% indSources = [indCygDec, indCygRa ; indCasDec, indCasRa];
%---------------------------------------------------------------------%


for n = 1:nSources
    % [maxValue,indMax] = max(abs(intensityClean(:)));
    % [indMaxDec, indMaxRa] = ind2sub(size(intensityMF),indMax);
    maxValue = intensityClean(indSources(n,1),indSources(n,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% When 2DFFT Imaging is used
    %---------------------------------------------------------------------%
    PSF = calculatePointSpreadFunction2DFFT(tgExact,indSources(n,1),indSources(n,2),504,1,1024,447);
    %---------------------------------------------------------------------%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% When Standard Imaging is used
    %---------------------------------------------------------------------%
%     PSF = calculatePointSpreadFunction(tgExact,indSources(n,1),indSources(n,2),1,1024);
    %---------------------------------------------------------------------%
    
    PSF = PSF * maxValue/max(abs(PSF(:)));
    
    intensityClean = intensityClean - PSF;
end

mesh(rA,dec,real(intensityDirty),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);

figure
mesh(rA,dec,real(intensityClean),'FaceColor','interp','LineStyle','none')
xlabel('RA (Hour)')
ylabel('Dec (Degrees)')
view(2)
colormap('jet');
set(gca,'xlim',[min(rA) max(rA)]);