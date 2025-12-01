clear; close all; clc;
%% ------------------ USER INPUT (leave TLE as provided) ------------------
tle1 = '1 65053U 25163A   25332.14766872  .00000367  00000-0  12401-3 0  9999';
tle2 = '2 65053  98.4058 157.6931 0001225  89.6258 270.5072 14.42501166 17400';
% Ground station (change to your location if needed)
gs_lat = 13.02;    % degrees
gs_lon = 77.57;    % degrees
gs_alt_m = 900;    % meters
% Carrier frequency for Doppler calculation
fc = 437.5e6;      % Hz
% Time window
start_utc = datetime('now','TimeZone','UTC');   % or set explicitly
duration_s = 900;    % seconds
dt = 0.5;            % seconds
%% ----------------------------------------------------------------------
% Basic checks
if ~exist('twoline2rv','file')
    error('twoline2rv not found. Ensure SGP4 .m files are in current folder.');
end
if ~exist('sgp4','file')
    error('sgp4 not found. Ensure SGP4 .m files are in current folder.');
end
% Build times
times = start_utc + seconds(0:dt:duration_s);
M = numel(times);

% ------------------- DEFINITIVE twoline2rv BLOCK (FIXED) ----------------
satrec = [];
n_in = nargin('twoline2rv');
fprintf('twoline2rv reports %d input(s). Attempting initialization...\n', n_in);

try
    if n_in == 2
        % twoline2rv(tle1, tle2)
        satrec = twoline2rv(tle1, tle2);
        fprintf('twoline2rv called as (tle1, tle2)\n');

    elseif n_in >= 4
        % **FIXED CALL ORDER based on user's twoline2rv.m function header**
        % Expected signature: twoline2rv(longstr1, longstr2, opsmode, whichconst)
        opsmode = 'i';   % 'i' for Improved (or 'a' for AFSPC if 'i' fails)
        whichconst = 72; % WGS-72 constants

        satrec = twoline2rv(tle1, tle2, opsmode, whichconst);
        fprintf('twoline2rv called as (tle1, tle2, opsmode, whichconst)\n');
        
    else
        error('twoline2rv requires 2 or 4 input arguments, but %d were detected. Cannot proceed.', n_in);
    end

    % Final check after successful call attempt
    if isstruct(satrec) && isfield(satrec, 'error') && satrec.error ~= 0
        error('SGP4 Initialization Error (code %d): Check TLE syntax or SGP4 constants.', satrec.error);
    end

catch MEt
    error('twoline2rv attempt failed: %s', MEt.message);
end
% ------------------------------------------------------------------------

% capture epoch
jd_epoch = [];
if isfield(satrec,'jdsatepoch')
    jd_epoch = satrec.jdsatepoch;
    if isfield(satrec,'jdsatepochf')
        jd_epoch = jd_epoch + satrec.jdsatepochf;
    end
end

% Ground station ECEF (km) and velocity
[gs_x, gs_y, gs_z] = geodetic2ecef(gs_lat, gs_lon, gs_alt_m/1000);
gs_pos_ecef = [gs_x; gs_y; gs_z];
omega_earth = 7.2921150e-5; % rad/s
omega_vec = [0;0;omega_earth];
gs_vel_ecef = cross(omega_vec, gs_pos_ecef); % km/s

% Preallocate outputs
unix_ts = zeros(M,1);
fd_hz = zeros(M,1);
range_km = zeros(M,1);
c_km_s = 299792.458;

% Detect sgp4 nargout and start loop
n_out_sgp4 = nargout('sgp4');
fprintf('sgp4 reports %d output(s)\n', n_out_sgp4);
fprintf('Propagating %d epochs (dt=%.3fs)...\n', M, dt);

for i = 1:M
    t = times(i);
    unix_ts(i) = posixtime(t);

    % compute tsince (minutes)
    if ~isempty(jd_epoch)
        jd_now = juliandate(t);
        tsince_min = (jd_now - jd_epoch) * 1440;
    else
        % Fallback calculation if jd_epoch is not available in satrec
        if isfield(satrec,'epochyr') && isfield(satrec,'epochdays')
            yr = satrec.epochyr;
            if yr < 57, yr = yr + 2000; else yr = yr + 1900; end
            epoch_datetime = datetime(yr,1,0,'TimeZone','UTC') + days(satrec.epochdays);
            tsince_min = minutes(t - epoch_datetime);
        else
            tsince_min = 0; % Should not happen with valid TLE/satrec
        end
    end

    % Call sgp4
    r_eci = []; v_eci = [];
    try
        % Simplified robust sgp4 call based on common signatures
        if n_out_sgp4 >= 3
            % Assumes error code is first output: [error, pos, vel] = sgp4(...)
            [~, r_eci_temp, v_eci_temp] = sgp4(satrec, tsince_min); 
            r_eci = r_eci_temp(:); v_eci = v_eci_temp(:);
        else
            % Attempt the 2-output call (most common): [pos, vel] = sgp4(...)
            [r_eci_temp, v_eci_temp] = sgp4(satrec, tsince_min);
            r_eci = r_eci_temp(:); v_eci = v_eci_temp(:);
        end
    catch sgpe
        error(['sgp4 failed propagation at epoch %d. Error: %s'], i, sgpe.message);
    end

    % ensure column vectors and non-empty
    if isempty(r_eci) || isempty(v_eci) || numel(r_eci) < 3
        warning('SGP4 failed to compute position/velocity at time %s. Skipping point.', char(t));
        range_km(i) = NaN; fd_hz(i) = NaN;
        continue;
    end
    r_eci = r_eci(1:3); v_eci = v_eci(1:3);

    % Convert TEME/ECI -> ECEF approx via GST rotation
    gst = gstime_from_datetime(t);
    Rz = [cos(gst) sin(gst) 0; -sin(gst) cos(gst) 0; 0 0 1];
    sat_pos_ecef = Rz * r_eci;
    sat_vel_ecef = Rz * v_eci + cross(omega_vec, sat_pos_ecef);

    % Range & Range-Rate
    r_vec = gs_pos_ecef - sat_pos_ecef;
    r = norm(r_vec);
    range_km(i) = r;
    rel_vel = gs_vel_ecef - sat_vel_ecef;
    range_rate = dot(rel_vel, r_vec) / r; % km/s

    % Doppler (Hz)
    fd_hz(i) = - (range_rate / c_km_s) * fc;
end

% Save CSV
% FIX: Use times(:) to ensure the datetime array is a column vector
tbl = table(times(:), unix_ts, range_km, fd_hz, ...
    'VariableNames',{'utc_iso','unix_t','range_km','fd_hz'});
writetable(tbl,'doppler_compensation.csv');
fprintf('Wrote doppler_compensation.csv with %d rows\n', M);

% Quick Plot
figure('Name','Predicted Doppler from TLE');
plot((0:M-1)*dt, fd_hz,'LineWidth',1.4);
xlabel('Time (s)'); ylabel('Doppler (Hz)'); grid on;
title('Predicted Doppler vs time (from TLE via SGP4)');

%% ---------------- Helper functions (REQUIRES MATLAB TOOLBOXES) ----------------
% Note: The function juliandate used inside gstime_from_datetime requires
% the Mapping Toolbox or other date functions.
function [x,y,z] = geodetic2ecef(lat_deg, lon_deg, alt_km)
    a = 6378.137; e2 = 6.69437999014e-3;
    lat = deg2rad(lat_deg); lon = deg2rad(lon_deg);
    N = a ./ sqrt(1 - e2 .* sin(lat).^2);
    x = (N + alt_km) .* cos(lat) .* cos(lon);
    y = (N + alt_km) .* cos(lat) .* sin(lon);
    z = ((1 - e2) .* N + alt_km) .* sin(lat);
end

function gst = gstime_from_datetime(dt)
    % This function calculates Greenwich Sidereal Time (GST) in radians
    % from a datetime object (dt).
    % Note: Requires the 'juliandate' function, typically from the Mapping or Aerospace Toolbox.
    jd = juliandate(dt);
    T = (jd - 2451545.0)/36525;
    GMST_sec = 67310.54841 + (876600*3600 + 8640184.812866)*T + 0.093104*T^2 - 6.2e-6*T^3;
    GMST_sec = mod(GMST_sec,86400);
    if GMST_sec < 0, GMST_sec = GMST_sec + 86400; end
    gst = GMST_sec * (2*pi/86400);
end