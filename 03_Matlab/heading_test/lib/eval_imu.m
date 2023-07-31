clc, clear variables
%%

data = readmatrix('putty_09.log');

time = data(:,10) * 1e-3;
time = time - time(1);
Ts = median(diff(time));
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
ind_quat = 11:14;
ind_rpy  = 15:17;

figure(1)
plot(diff(time))

figure(2)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), time, data(:,ind_gyro) * 180/pi), grid on, ylabel('Gyro (deg/sec)')
ax(2) = nexttile;
plot(ax(2), time, data(:,ind_acc)), grid on, ylabel('Acc (m/s^2)')
ax(3) = nexttile;
plot(ax(3), time, data(:,ind_mag)), grid on, ylabel('Mag'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

Teval = [[ 2 3]; ...
         [ 5 7]; ...
         [ 9 10]];

figure(3)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), time, cumtrapz(time, data(:,ind_gyro)) * 180/pi), grid on, ylabel('Gyro Integral (deg)')
ax(2) = nexttile;
plot3(data(:,ind_acc(1)), data(:,ind_acc(2)), data(:,ind_acc(3)), 'k'), grid on, hold on
ind = time > Teval(1,1) & time < Teval(1,2);
plot3(data(ind,ind_acc(1)), data(ind,ind_acc(2)), data(ind,ind_acc(3)), 'b')
ind = time > Teval(2,1) & time < Teval(2,2);
plot3(data(ind,ind_acc(1)), data(ind,ind_acc(2)), data(ind,ind_acc(3)), 'color', [0 0.5 0])
ind = time > Teval(3,1) & time < Teval(3,2);
plot3(data(ind,ind_acc(1)), data(ind,ind_acc(2)), data(ind,ind_acc(3)), 'r'), hold off
axis([-1 1 -1 1 -1 1] * 10), axis equal
xlabel('x-Axis'), ylabel('y-Axis'), zlabel('z-Axis')
ax(3) = nexttile;
plot3(data(:,ind_mag(1)), data(:,ind_mag(2)), data(:,ind_mag(3)), 'k'), grid on, hold on
ind = time > Teval(1,1) & time < Teval(1,2);
plot3(data(ind,ind_mag(1)), data(ind,ind_mag(2)), data(ind,ind_mag(3)), 'b')
ind = time > Teval(2,1) & time < Teval(2,2);
plot3(data(ind,ind_mag(1)), data(ind,ind_mag(2)), data(ind,ind_mag(3)), 'color', [0 0.5 0])
ind = time > Teval(3,1) & time < Teval(3,2);
plot3(data(ind,ind_mag(1)), data(ind,ind_mag(2)), data(ind,ind_mag(3)), 'r'), hold off
axis([-1 1 -1 1 -1 1] * 0.5), axis equal

data_spectras = data(:,[ind_gyro, ind_acc, ind_mag]);
data_spectras = data_spectras - mean(data_spectras);
[pxx, f] = pwelch(data_spectras, round(1/Ts), [], [], 1/Ts, 'Power');
spectras = sqrt(pxx*2); % power -> amplitude (dc needs to be scaled differently)

figure(4)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), f, spectras(:,ind_gyro)), grid on, set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylabel('Gyro (deg/sec)')
ax(2) = nexttile;
plot(ax(2), f, spectras(:,ind_acc)), grid on, set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylabel('Acc (m/s^2)')
ax(3) = nexttile;
plot(ax(3), f, spectras(:,ind_mag)), grid on, set(gca, 'XScale', 'log'), set(gca, 'YScale', 'log'), ylabel('Mag'), xlabel('Frequency (Hz)')
linkaxes(ax, 'x'), clear ax, xlim([0 1/2/Ts])

figure(5)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(1)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(data(:,ind_rpy(1)))  * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(2) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(2)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(data(:,ind_rpy(2)))  * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(3) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(3)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(data(:,ind_rpy(3)))  * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])
