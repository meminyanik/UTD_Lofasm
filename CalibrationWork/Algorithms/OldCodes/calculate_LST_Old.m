function [GST_HMS,GST1_HMS,GST2_HMS,GST4_HMS,LST_HMS] = calculate_LST_Old(longitude,hour,minute,second,day,month,year,isLeapYear)

% GSTx_HMS: Greenwich Siderial Time (GST) in hours, minutes, seconds format
% LST_HMS: Local Siderial Time (LST) in hours, minutes, seconds format

% Example Function Call
% Find the number of days from J2000.0 for 05:45:11 hrs UT on 2016 July 21th
% [GST_HMS,GST1_HMS,GST2_HMS,GST4_HMS,LST_HMS] = calculate_LST(-116.791509,5,45,11,21,7,2016,leapyear(2016))

 % for debugging
 % hour = 5; minute = 45; second = 11; day = 21; month = 7;
 % year = 2016; isLeapYear = leapyear(2016);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LoFASM IV Antenna Coordinates
%-------------------------------------------------------------------------%
% mainLong1Map = -116.793287;
% mainLong2Map = -116.793316;
% mainLong3Map = -116.793300;
% mainLong4Map = -116.793251;
% mainLong5Map = -116.793221;
% mainLong6Map = -116.793239;
%
% outLongMap = -116.791509;
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the days from J2000
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
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Julian Date
% jd = juliandate([y,mo,d,h,mi,s])
% JD = juliandate([year,month,day,hour,minute,second]);
% J2000 refers to the instant of January 1, 2000, 11:58:55.816 UTC
% numOfDays = juliandate(year,month,day,hour,minute,second) - juliandate(2000,1,1,11,58,55.816);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Find the Greenwich Siderial Time (GST) for Old Method

% Find the fraction hours
startHour = hour + minute/60 + second/3600;

% Calculate J2000 days
dayNum = day;

if isLeapYear
    daysMonth = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
else
    daysMonth = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
end
% daysYear = [1998, 1999, ..., 2021]
daysYear = [-731.5, -366.5, -1.5, 364.5, 729.5, 1094.5, 1459.5, 1825.5, 2190.5, 2555.5, 2920.5, 3286.5, ...
    3651.5, 4016.5, 4381.5, 4747.5, 5112.5, 5477.5, 5842.5, 6208.5, 6573.5, 6938.5, 7303.5, 7669.5];

daysSinceJanuary = daysMonth(month);
daysSinceJ2000 = daysYear(year - 1997);
numOfDays = daysSinceJ2000 + daysSinceJanuary + dayNum;

GST1 = 100.46 + 0.985647 * numOfDays + 15*startHour;
GST1 = mod(GST1,360)/15;
GST1_HMS = decimalHours2HoursMinutesSeconds(GST1);
%-------------------------------------------------------------------------%

JD2000 = juliandate(2000,1,1,12,0,0); % Or juliandate(2000,1,1,11,58,55.816);
% JD2000 = 2451545.0;

%-------------------------------------------------------------------------%
% Find the Greenwich Siderial Time (GST) for Second Method
numOfDays = juliandate(year,month,day) - JD2000;
GST2 = 100.46 + 0.985647 * numOfDays + 15*startHour;
GST2 = mod(GST2,360)/15;
GST2_HMS = decimalHours2HoursMinutesSeconds(GST2);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Find the Greenwich Siderial Time (GST) for Third Method
JD = juliandate(year,month,day);
UT = hour + minute/60 + second/3600;
T = (JD - JD2000)/36525.0;
GST = 6.697374558 + (2400.051336*T)+(0.000025862*T^2)+(UT*1.0027379093);
GST = mod(GST,24);
GST_HMS = decimalHours2HoursMinutesSeconds(GST);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Find the Greenwich Siderial Time (GST) for Fourth Method
JD = juliandate(year,month,day,hour,minute,second);
GST4 = 18.697374558 + 24.06570982441908 * (JD - JD2000);
GST4 = mod(GST4,24);
GST4_HMS = decimalHours2HoursMinutesSeconds(GST4);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Find the Local Siderial Time (LST) of Antenna
LST = mod(GST + longitude/15,24);
LST_HMS = decimalHours2HoursMinutesSeconds(LST);
%-------------------------------------------------------------------------%

    function hms = decimalHours2HoursMinutesSeconds(h)
        HH = floor(h);
        MM = floor((h-HH)*60);
        SS = (h-HH)*3600 - MM*60;
        hms = HH+"h " + MM+"m " + SS+"s";
    end

end

