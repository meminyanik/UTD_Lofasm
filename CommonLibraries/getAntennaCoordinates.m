% This code block calculates antenna baseline coordinates of staaions

% Before running this code, following parameters should be defined
% - stationID: Indicates LoFASM Station ID

% Output: Antenna coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Antenna Coordinates
%-------------------------------------------------------------------------%
if stationID == 1
    inMainLat2Map = 26.555430;
    inMainLong2Map = -97.442018;
    inMainLat5Map = 26.555451;
    inMainLong5Map = -97.441933;
    
    outMainLat2Map = 26.555394;
    outMainLong2Map = -97.442031;
    outMainLat5Map = 26.555489;
    outMainLong5Map = -97.441920;
    
    mainAltMap = 3;
    
    outLatMap = 26.556041;
    outLongMap = -97.440489;
    outAltMap = 3;
    
elseif stationID == 4
    % mainLat1Map = 35.247190;
    % mainLong1Map = -116.793287;
    inMainLat2Map = 35.247222;
    inMainLong2Map = -116.793316;
    % mainLat3Map = 35.247256;
    % mainLong3Map = -116.793300;
    % mainLat4Map = 35.247264;
    % mainLong4Map = -116.793251;
    inMainLat5Map = 35.247233;
    inMainLong5Map = -116.793221;
    % mainLat6Map = 35.247196;
    % mainLong6Map = -116.793239;
    mainAltMap = distdim(3547,'ft','m');
    
    outLatMap = 35.247469;
    outLongMap = -116.791509;
    uDistance = 3; % Zenith baseline in meters
    outAltMap = distdim(3525,'ft','m') + uDistance;
end
%-------------------------------------------------------------------------%