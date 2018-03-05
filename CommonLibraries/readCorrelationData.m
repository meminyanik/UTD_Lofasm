%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
currentPath = pwd();
addpath(genpath([currentPath,'/../../../Algorithms']));
addpath(genpath([currentPath,'/../../../RecordedData']));

% This code block reads recorded correlation data from harddisk

% Before running this code, following parameters should be defined
% - Simulation Mode: flag indicates simulation mode
% - nSAv: number of Samples for averaging
% - stationID: Indicates LoFASM Station ID
% - Whitening: flag indicates whitening

% Output: RData

%--------------------------------------------------------------------------
% Load Data
%-------------------------------------------------------------------------%
if SimulationMode
    %--------------------------------------------------------------------------
    % Load Simulated Data
    %-------------------------------------------------------------------------%
    RlTerm = ABSim;
    
else
    %--------------------------------------------------------------------------
    % Load Real Data
    %-------------------------------------------------------------------------%
    if nSAv == 447
        if (stationID == 1)
            load('L1_AClTerm');
            RlTerm = L1_AClTerm.Data;
        elseif (stationID == 4)
            load('ABlTermPZ');
            RlTerm = ABlTerm;
        end
        if Whitening
            RlTerm = whiteningApproach1(RlTerm);
        end
    elseif nSAv == 200
        if (stationID == 1)
            load('L1_AClTerm_w'); % Data with Whitening first
            RlTerm = L1_AClTerm_w.Data;
        elseif (stationID == 4)
            % TBD
        end
    else
        if (stationID == 1)
            load('L1_AC.mat')
            R = L1_AC;
        elseif (stationID == 4)
            % TBD
        end
        
        if Whitening
            R = whiteningApproach1(R);
        end
        
        RlTerm = longTermAverage(R,nSAv);
    end
end

RlT = RlTerm; %./ (AAlTerm+BBlTerm); % clear ABlTerm