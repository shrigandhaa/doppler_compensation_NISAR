%--------------------------------------------------------------------------
%                   SGP4 Orbit Propagator
%
% References:
% Hoots, Felix R., and Ronald L. Roehrich. 1980. Models for Propagation of
% NORAD Element Sets. Spacetrack Report #3. U.S. Air Force: Aerospace Defense
% Command.
% 
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
% , New York; 4th edition (2013).
% 
% Last modified:   2023/08/30   Meysam Mahooti
%--------------------------------------------------------------------------
clc
clear all
format long g

global const % Astronomical Constants
SAT_Const

opsmode= 'a'; % afspc approach
whichconst = 72;


% TLE file name
fname = 'GRACE-FO1.txt';

% Open the TLE file and read it
fid = fopen(fname, 'r');

sat_name = fgetl(fid); % read satellite's name
longstr1 = fgetl(fid); % read first line
longstr2 = fgetl(fid); % read second line
satrec = twoline2rv(longstr1, longstr2, opsmode, whichconst);

tsince = 1440; % amount of time in which you are going to propagate satellite's state vector forward (+) or backward (-) [minutes] 
% call the propagator to get the final state vector value
[satrec,rteme,vteme] = sgp4(satrec, tsince);

fprintf('Ephemeris in the True Equator Mean Equinox coordinate system based on the epoch of the specified TLE.\n');
fprintf('     TSINCE              X                Y                Z     [km]\n');
fprintf(' %9.1f%22.8f%18.8f%18.8f \n',tsince,rteme(1),rteme(2),rteme(3));
fprintf('                       XDOT             YDOT             ZDOT    [km/s]\n');
fprintf('  %28.8f%18.8f%18.8f \n\n',vteme(1),vteme(2),vteme(3));

% read Earth orientation parameters
fid = fopen('EOP-All.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
while ~feof(fid)
    tline = fgetl(fid);
    k = strfind(tline,'NUM_OBSERVED_POINTS');
    if (k == 1)
        numrecsobs = str2num(tline(21:end));
        tline = fgetl(fid);
        for i=1:numrecsobs
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        for i=1:4
            tline = fgetl(fid);
        end
        numrecspred = str2num(tline(22:end));
        tline = fgetl(fid);
        for i=numrecsobs+1:numrecsobs+numrecspred
            eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
        end
        break
    end
end
fclose(fid);

MJD_UTC = (satrec.jdsatepoch+satrec.jdsatepochf-2400000.5)+tsince/1440;

% Earth Orientation Parameters
[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
MJD_UT1 = MJD_UTC + UT1_UTC/86400;
MJD_TT  = MJD_UTC + TT_UTC/86400;
T = (MJD_TT-const.MJD_J2000)/36525;
[reci, veci] = teme2eci(rteme',vteme',T,dpsi,deps)
[recef,vecef] = teme2ecef(rteme',vteme',T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2)
[rtod, vtod] = ecef2tod(recef,vecef,T,MJD_UT1+2400000.5,LOD,x_pole,y_pole,2,dpsi,deps)

