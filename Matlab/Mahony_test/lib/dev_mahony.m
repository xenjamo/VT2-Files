clc, clear variables
%%

data = readmatrix('putty_14.log'); % gyro and acc
% data = readmatrix('putty_15.log'); % gyro, acc and mag

time = data(:,10) * 1e-3;
time = time - time(1);
Ts = median(diff(time));
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
ind_quat = 11:14;
ind_rpy  = 15:17;

figure(1)
plot(diff(time)), grid on

figure(2)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), time, data(:,ind_gyro) * 180/pi), grid on, ylabel('Gyro (deg/sec)')
ax(2) = nexttile;
plot(ax(2), time, data(:,ind_acc)), grid on, ylabel('Acc (m/s^2)')
ax(3) = nexttile;
plot(ax(3), time, data(:,ind_mag)), grid on, ylabel('Mag'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

figure(3)
plot(time, cumtrapz(time, data(:,ind_gyro)) * 180/pi), grid on
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
xlim([0, max(time)])

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

figure(4)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
ang_ind = 1;
plot(time, angleFun(cumtrapz(time, data(:,ind_gyro(ang_ind)))) * 180/pi), grid on, hold on
plot(time, angleFun(rpyRP (:,ang_ind)) * 180/pi)
plot(time, angleFun(rpyRPY(:,ang_ind)) * 180/pi)
plot(time, angleFun(data(:,ind_rpy(ang_ind))) * 180/pi), hold off
ylabel('Roll (deg)'), xlabel('Time (sec)')
ax(2) = nexttile;
ang_ind = 2;
plot(time, angleFun(cumtrapz(time, data(:,ind_gyro(ang_ind)))) * 180/pi), grid on, hold on
plot(time, angleFun(rpyRP (:,ang_ind)) * 180/pi)
plot(time, angleFun(rpyRPY(:,ang_ind)) * 180/pi)
plot(time, angleFun(data(:,ind_rpy(ang_ind))) * 180/pi), hold off
ylabel('Pitch (deg)'), xlabel('Time (sec)')
ax(3) = nexttile;
ang_ind = 3;
plot(time, angleFun(cumtrapz(time, data(:,ind_gyro(ang_ind)))) * 180/pi), grid on, hold on
plot(time, angleFun(rpyRP (:,ang_ind)) * 180/pi)
plot(time, angleFun(rpyRPY(:,ang_ind)) * 180/pi)
plot(time, angleFun(data(:,ind_rpy(ang_ind))) * 180/pi), hold off
ylabel('Yaw (deg)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

figure(5)
tiledlayout_ = tiledlayout(2,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(time, biasRP  * 180/pi), grid on
ylabel('Gyro Bias (deg/s)')
ax(2) = nexttile;
plot(time, biasRPY * 180/pi), grid on
ylabel('Gyro Bias (deg/s)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

% figure(6)
% subplot(211)
% plot(time, data(:,ind_quat)), grid on
% ylabel('Quaternion')
% subplot(212)
% plot(time, quat_), grid on
% ylabel('Quaternion'), xlabel('Time (sec)')
