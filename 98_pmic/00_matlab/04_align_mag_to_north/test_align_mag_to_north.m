clc, clear variables
addpath ..\99_fcn_bib\
%% 
load data.mat

data = [gyro, accel, mag_calib, t, quat, rpy, headMot];

time = data(:, 10);
time = time - time(1);
Ts = median( diff(time) );
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
ind_quat = 11:14;
ind_rpy  = 15:17;
ind_headMot = 18;

%%

% Sources to get magnetic declination and inclination:

% - https://earth.google.com/web/search/Winterthur/@47.48925878,8.74035197,446.306308a,242.88835455d,35y,0h,0t,0r/data=CigiJgokCYI4KOL8OzhAEXw4KOL8OzjAGSsiU3ylgkVAITUWFRALAFDA

% - https://www.magnetic-declination.com/
% 02.06.2023:
% You clicked here:
% Latitude: 47° 29' 21.1" N
% Longitude: 8° 44' 25.7" E
% WINTERTHUR
% Magnetic Declination: +3° 22'
% Declination is POSITIVE (EAST)
% Inclination: 63° 34'
% Magnetic field strength: 48314.4 nT

% - https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#declination
% 02.06.2023:
% 2023-06-02 	3° 21' E  ± 0° 22'  changing by  0° 10' E per year

% Info:
% - https://www.ngdc.noaa.gov/geomag/geomaginfo.shtml

% short rotation
myWrapToPi = @(ang) atan2( sin(ang), cos(ang) );

Tmean = 2.0; % make sure we have proper data here
ind_avg = time < Tmean;

yaw_mag  = -atan2( -data(:,ind_mag(2) ), ...
                  -data(:,ind_mag(1) ) );
yaw_mag0 = median(yaw_mag(ind_avg));
fprintf("   yaw_mag0 (deg): %0.2f\n", yaw_mag0 * 180/pi);

mag_decli_online = dms2degrees( [ 3, 22, 0] ) * pi/180;
mag_incli_online = dms2degrees( [63, 34, 0] ) * pi/180;
fprintf("   mag_decli_online (deg): %0.2f\n", mag_decli_online * 180/pi);
fprintf("   mag_incli_online (deg): %0.2f\n", mag_incli_online * 180/pi);

mag0 = mean( data(ind_avg, ind_mag) ).';
mag0_xy_proj = mag0; mag0_xy_proj(3) = 0;
mag_incli_meas = calcAngleBetween2Vectors(mag0, mag0_xy_proj);
fprintf("   mag_inclination_meas (deg): %0.2f\n", mag_incli_meas * 180/pi);

yaw_mag = yaw_mag - mag_decli_online;

figure(1)
plot(time, myWrapToPi( [yaw_mag, -data(:,ind_headMot)] ) * 180/pi),  grid on
xlabel('Time (sec)'), ylabel('Yaw (deg)')
legend('Magnetometer (this is an approximation)', ...
       'GNSS headMot (assume we are purely flying forward)', ...
       'location', 'best')
xlim([0 100])

figure(2)
arrow_length = 1.0;
quiver3(0, 0, 0, arrow_length, 0, 0, 'LineWidth', 2, 'AutoScale', 'off', 'color', [1 0 0]), hold on
quiver3(0, 0, 0, 0, arrow_length, 0, 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0.5 0])
quiver3(0, 0, 0, 0, 0, arrow_length, 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0 1])
quiver3(0, 0, 0, mag0(1), mag0(2), mag0(3), 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 1 1]), hold off
xlabel('x-Axis'), ylabel('y-Axis'), zlabel('z-Axis'), axis equal
view(-133, 50)

