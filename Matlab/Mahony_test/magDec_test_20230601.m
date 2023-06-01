clc, clear variables

%% 
load data.mat

data = [gyro, accel, mag_calib, t, quat, rpy, headMot];

time = data(:,10);
time = time - time(1);
Ts = median(diff(time));
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
ind_quat = 11:14;
ind_rpy  = 15:17;
ind_headMot = 18;

%% 
Tmean = 1.0; % ???
ind_avg = time < Tmean;
mag0 = mean(data(ind_avg,ind_mag)).';

yaw_mag0 = median( atan2(data(ind_avg,ind_mag(2)), data(ind_avg,ind_mag(1))) );

yaw_mag = atan2(data(:,ind_mag(2)), data(:,ind_mag(1)));
yaw_mag = wrapTo2Pi(yaw_mag - 0*yaw_mag0);

figure(839)
plot(time, [yaw_mag,data(:,ind_headMot)] * 180/pi)
grid on;

% R = Rz(psi0)^T * Rz(psi)*Ry(theta)*Rx(phi)