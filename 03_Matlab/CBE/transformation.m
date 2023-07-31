clear all; clc

load data.mat;

data = [gyro, accel, mag_calib, t, quat, rpy];
[m,n] = size(t);

time = data(:,10);
time = time - time(1);
Ts = median(diff(time));
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;

% bessel
p = 2;         % pole at p rad/s
kp = 2 * p;
ki = kp^2 / 3;

para.kp = kp;
para.ki = ki;

rpy0 = 0 * [60, -60, 0] * 180/pi;
quat0 = rpy2quat(rpy0).';

[quatRP , biasRP ] = mahonyRP (data(:,ind_gyro), data(:,ind_acc), para, Ts, quat0);
[quatRPY, biasRPY] = mahonyRPY(data(:,ind_gyro), data(:,ind_acc), data(:,ind_mag), para, Ts, quat0);

rpyRP  = quat2rpy(quatRP );
rpyRPY = quat2rpy(quatRPY);

angleFun = @wrapToPi;

%% transformation matrix

CBE = rotmat(quaternion(quatRPY),'frame');
CEB = rotmat(quaternion(quatRPY),'point');
accelBE = zeros(m,3);

for i = 1:m
    accelBE(i,:) = accel(i,:)*CBE(:,:,i);
end

figure(1)
subplot(211)
plot(time, accel)
grid on;
ylabel("accel body frame")
subplot(212)
plot(time, accelBE)
grid on;
ylabel("accel world frame")