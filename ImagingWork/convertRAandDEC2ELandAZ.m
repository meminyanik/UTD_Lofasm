%% Intro
% http://www.stargazing.net/kepler/altaz.html#twig04
% You will find out how to calculate the Azimuth (AZ) and Elevation (EL) of an object in the sky
% if you know the date, time (UT) and the location of your observing site
% together with the Right Ascension (RA) and Declination (DEC) of the object

% As a concrete example, calculate the EL and AZ of the Messier object M13 for 10th August 1998 at 2310 hrs UT for Birmingham UK

% The RA and DEC of M13 are given as;
%  RA  = 16 h 41.7 min
%  DEC = 36 d 28   min

% Latitude and longitude of Birmingham UK as;
%  LAT = 52 d 30 min North
%  LONG = 1 d 55 min West

% We will need these figures in decimal form, along with the time;
%  RA   = 16 h 41.7 min     = 16 + 41.7/60 = 16.695     hrs
%  DEC  = 36 d 28   min     = 36 + 28/60   = 36.466667  degs
%  Time = 2310 hrs          = 23 + 10/60   = 23.166667  hrs
%  LAT  = 52 d 30 min North = 52 + 30/60   = 52.5       degs
%  LONG =  1 d 55 min West  = -(1 + 55/60) = -1.9166667 degs

% Longitudes west are counted as negative, and East counted as positive
% RA  = 16.695 * 15 = 250.425 degrees

% Calculate the Local Siderial Time (LST), and then work out the Hour Angle (HA) of the object
% Then we can use some standard formulas from spherical trigonometry to transform the HA and DEC to the ALT and AZ.
%%

%% Calculating the days from J2000
% -------------------------------
% 
% The tables below can be used to calculate the number of days and
% the fraction of a day since the epoch J2000. If you need the
% number of Julian centuaries, then just divide the 'day number' 
% by 36525.
% 
% Table A                |  Table B
% Days to beginning of   |  Days since J2000 to
% month                  |  beginning of each year
%                        |
% Month   Normal   Leap  |  Year   Days    |  Year   Days
%         year     year  |                 |
%                        |                 |
% Jan       0        0   |  1998   -731.5  |   2010  3651.5
% Feb      31       31   |  1999   -366.5  |   2011  4016.5
% Mar      59       60   |  2000     -1.5  |   2012  4381.5
% Apr      90       91   |  2001    364.5  |   2013  4747.5
% May     120      121   |  2002    729.5  |   2014  5112.5    
% Jun     151      152   |  2003   1094.5  |   2015  5477.5
% Jul     181      182   |  2004   1459.5  |   2016  5842.5
% Aug     212      213   |  2005   1825.5  |   2017  6208.5
% Sep     243      244   |  2006   2190.5  |   2018  6573.5
% Oct     273      274   |  2007   2555.5  |   2019  6938.5
% Nov     304      305   |  2008   2920.5  |   2020  7303.5
% Dec     334      335   |  2009   3286.5  |   2021  7669.5
% 
% Worked Example
%                                                                
% To find the number of days from J2000.0 for 2310 hrs UT on
% 1998 August 10th, do the following;
% 
% 1. divide the number of minutes by 60 to obtain the decimal
% fraction of an hour, here 10/60 = 0.1666667
%    
% 2. add this to the hours, then divide the total by 24 to obtain
% the decimal fraction of the day, here 23.166667/24 = 0.9652778 
% This is the first number used below
%    
% 3. find from table A the number of days to the beginning of
% August from the start of the year, here 212 days
%    
% 4. write down the day number within the month, here 10 above
% 
% 5. find from table B the days since J2000.0 to the beginning of
% the year, here -731.5
%    
% 6. add these four numbers.
% 
% For the date above;
% 
%    0.9652778 + 212 + 10 - 731.5 = -508.53472 days from J2000.0
% 
% Note that dates which fall before J2000.0 will have negative day
% numbers. Keep the negative sign in any calculations. 
%%
% Calculate Julian date for 2016 July 21th, at 05:45:11 a.m.:
% mjd = mjuliandate(y,mo,d,h,mi,s)
% jd = juliandate(y,mo,d,h,mi,s)
% since noon Universal Time on January 1, 4713 BCE
% jd = juliandate(2016,07,21,05,45,11);

% Find the number of days from J2000.0 for 05:45:11 hrs UT on 2016 July 21th,
startHour = 05 + 45/60 + 11/3600;
% startHour = startHour - 1;
fractionDay = startHour / 24;
dayNum = 21;
daysSinceJanuary = 182; %July
daysSinceJ2000 = 5842.5;
numOfDays = daysSinceJ2000 + daysSinceJanuary + dayNum + fractionDay;

%% Local Siderial Time
% LST = 100.46 + 0.985647 * d + long + 15*UT
%   d    is the days from J2000, including the fraction of a day
%   UT   is the universal time in decimal hours
%   long is your longitude in decimal degrees, East positive      
% Add or subtract multiples of 360 to bring LST in range 0 to 360 degrees.

% Find the local siderial time for 2310 UT, 10th August 1998
% at Birmingham UK (longitude 1 degree 55 minutes west).
% 
% I know that UT = 23.166667
%              d = -508.53472 (last section)
%           long = -1.9166667  (West counts as negative)
% 
% so 
% 
% LST = 100.46 + 0.985647 * d + long + 15*UT
%     = 100.46 + 0.985647 * -508.53472 - 1.9166667 + 15 * 23.166667
%     = -55.192383 degrees
%     = 304.80762 degrees
%     
% note how we added 360 to LST to bring the number into the range
% 0 to 360 degrees.
%%
mainLat = 35 + 14/60 + 50.01/3600;
mainLong = -(116 + 47/60 + 35.78/3600);
mainLatMap = 35.247227;
mainLongMap = -116.793272;

%% myLong = -96.826593;

LST = 100.46 + 0.985647 * numOfDays + mainLongMap + 15*startHour;
LST = mod(LST,360);
% LST = wrapTo180(LST);

raCygnusA = (19 + 59/60 + 28.3/3600) * 15;
decCygnusA = 40 + 44/60 + 02/3600;
% Cygnus A: 19h59m28.3s +40d44m02s
raCasA = (23 + 23/60 + 27.9/3600) * 15;
decCasA = 58 + 48/60 + 42/3600;
% Cas A: 23h23m27.9s +58d48m42s

haCygnusA = LST - raCygnusA;
haCygnusA = mod(haCygnusA,360);
haCasA = LST - raCasA;
haCasA = mod(haCasA,360);

[elCygnusA,azCygnusA] = pos2elaz(mainLat,haCygnusA,decCygnusA);
% elCygnusA = elCygnusA*pi/180;
% azCygnusA = azCygnusA*pi/180;

[elCasA,azCasA] = pos2elaz(mainLat,haCasA,decCasA);
% elCasA = elCasA*pi/180;
% azCasA = azCasA*pi/180;

tSize = 504;
hRecording = (0*447*0.083886 : 447*0.083886 : (tSize-1)*447*0.083886)/3600;
LSTRecording = LST + 15*hRecording;
LSTRecording = mod(LSTRecording,360);

haCygnusARecording = LSTRecording - raCygnusA;
haCygnusARecording = mod(haCygnusARecording,360);
haCasARecording = LSTRecording - raCasA;
haCasARecording = mod(haCasARecording,360);

[elCygnusARecording,azCygnusARecording] = pos2elaz(mainLat,haCygnusARecording,decCygnusA);
[elCasARecording,azCasARecording] = pos2elaz(mainLat,haCasARecording,decCasA);

for n = 1:tSize
    if azCygnusARecording(n)>180
        azCygnusARecording(n) = azCygnusARecording(n)-360;
    end
    if azCasARecording(n)>180
        azCasARecording(n) = azCasARecording(n)-360;
    end
end

% elCygnusA = zeros(1,length(Time));
% azCygnusA = zeros(1,length(Time));
% elCasA = zeros(1,length(Time));
% azCasA = zeros(1,length(Time));
% 
% for n = 1:length(Time)
%     [elCygnusA(n),azCygnusA(n)] = pos2elaz(mainLat,Time(n),mainLong);
% end


function [el,az] = pos2elaz(latStation,hourAngleObject,decSource)

    % http://www.stargazing.net/kepler/altaz.html
    
    % If inputs are in degrees
    latStation = latStation*pi/180;
    hourAngleObject = hourAngleObject*pi/180;
    decSource = decSource*pi/180;

    % note: input args are in radians, and can be matrices
    el = asin(sin(decSource).*sin(latStation) + cos(decSource).*cos(latStation).*cos(hourAngleObject));
    az  = (sin(decSource) - sin(el).*sin(latStation)) ./ (cos(el).*cos(latStation));
    
    % remove roundoff errors
    [i,j] = find(abs(az)>1.0); 
    az(i,j) = sign(az(i,j)); 

    % change hemisphere
    az = acos(az);
    az(sin(hourAngleObject)>0) = deg2rad(360) - az(sin(hourAngleObject)>0); 
    
    el = rad2deg(el);
    az = rad2deg(az);
end