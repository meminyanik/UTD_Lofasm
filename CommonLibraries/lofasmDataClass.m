classdef lofasmDataClass
    properties
        stationID % 1,2,3,4
        arrayConfig % 0: inner ring - outtrigger, % 1: outer ring - outtrigger, % 2: inner ring - outer ring
        dataRecordTime % [h,m,s,d,m,y]
        nSAv % Number of samples used for averaging. nSAv=1 means no averaging
        whiteningApplied % 0: no whitening applied to data, 1: whitening applied before averaging
        Data; % Lofasm Data
    end
    
    methods
        % Object Constructor
        function lofasmData = lofasmDataClass(stationID,arrayConfig,dataRecordTime,nSAv,whiteningApplied,Data)
            if (nargin == 6)
                lofasmData.stationID = stationID;
                lofasmData.arrayConfig = arrayConfig;
                lofasmData.dataRecordTime = dataRecordTime;
                lofasmData.nSAv = nSAv;
                lofasmData.whiteningApplied = whiteningApplied;
                lofasmData.Data = Data;
            else
                error("Number of inputs is less than 6, please check again inputs")
            end
        end
    end
    
end