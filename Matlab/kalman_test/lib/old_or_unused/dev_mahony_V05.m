clc, clear variables
%%

% data = readmatrix('putty_09.log');
% data = readmatrix('putty_10.log');
data = readmatrix('putty_11.log');

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

figure(3)
plot(time, cumtrapz(time, data(:,ind_gyro)) * 180/pi), grid on
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')

% bessel
p = 2;         % pole at p rad/s
kp = 2 * p;
ki = kp^2 / 3;
N = size(data, 1);
quat = eye(4,1);
bias = zeros(3,1);
Eva = [0 0 1].';
quat_ = zeros(N, 4);
bias_ = zeros(N, 3);

% % betaflight (rp)
% for i = 1:N
%     g_n = [         2*(quat(2)*quat(4) - quat(1)*quat(3)); ...
%                     2*(quat(3)*quat(4) + quat(1)*quat(2)); ...
%             quat(1)^2 - quat(2)^2 - quat(3)^2 + quat(4)^2];
% 
%     acc_n = data(i, ind_acc).';
%     acc_n = acc_n ./ norm(acc_n);
%     e = cross(acc_n, g_n);
%     
%     bias = bias + ki * e * Ts;
% 
%     Q = [[-quat(2), -quat(3), -quat(4)]; ...
%          [ quat(1), -quat(4),  quat(3)]; ...
%          [ quat(4),  quat(1), -quat(2)]; ...
%          [-quat(3),  quat(2),  quat(1)]];
% 
%     dquat = Ts * 0.5 * Q * ( data(i,ind_gyro).' + bias + kp * e );
%     quat = quat + dquat;
%     quat = quat ./ norm(quat);
% 
%     quat_(i,:) = quat.';
%     bias_(i,:) = bias.';
% end

% betaflight (rpy)
for i = 1:N

    mag_n = data(i, ind_mag).';
    mag_n = mag_n ./ norm(mag_n);

%     // For magnetometer correction we make an assumption that magnetic field is perpendicular to gravity (ignore Z-component in EF).
%     // This way magnetic field will only affect heading and wont mess roll/pitch angles
% 
%     // (hx; hy; 0) - measured mag field vector in EF (assuming Z-component is zero)
%     // (bx; 0; 0) - reference mag field vector heading due North in EF (assuming Z-component is zero)
%     const float hx = rMat[0][0] * mx + rMat[0][1] * my + rMat[0][2] * mz;
%     const float hy = rMat[1][0] * mx + rMat[1][1] * my + rMat[1][2] * mz;
%     const float bx = sqrtf(hx * hx + hy * hy);

    rMat = quat2CBE(quat).';
    hx = rMat(1,:) * mag_n;
    hy = rMat(2,:) * mag_n;
    bx = sqrt(hx * hx + hy * hy);
    ez_ef = -(hy * bx);
    e = rMat(3,:).' * ez_ef;
    
    g_n = [         2*(quat(2)*quat(4) - quat(1)*quat(3)); ...
                    2*(quat(3)*quat(4) + quat(1)*quat(2)); ...
            quat(1)^2 - quat(2)^2 - quat(3)^2 + quat(4)^2];

    acc_n = data(i, ind_acc).';
    acc_n = acc_n ./ norm(acc_n);
    e = e + cross(acc_n, g_n);
    
    bias = bias + ki * e * Ts;

    Q = [[-quat(2), -quat(3), -quat(4)]; ...
         [ quat(1), -quat(4),  quat(3)]; ...
         [ quat(4),  quat(1), -quat(2)]; ...
         [-quat(3),  quat(2),  quat(1)]];

    dquat = Ts * 0.5 * Q * ( data(i,ind_gyro).' + bias + kp * e );
    quat = quat + dquat;
    quat = quat ./ norm(quat);

    quat_(i,:) = quat.';
    bias_(i,:) = bias.';
end

rpy_ = quat2rpy(quat_);

[z_, y_, x_] = quat2angle(quat_, "ZYX");
rpy_matlab = fliplr([z_, y_, x_]);

figure(4)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(1)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(rpy_(:,1)) * 180/pi)
plot(time, wrapToPi(data(:,ind_rpy(1))) * 180/pi)
plot(time, wrapToPi(rpy_matlab(:,1)) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(2) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(2)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(rpy_(:,2)) * 180/pi)
plot(time, wrapToPi(data(:,ind_rpy(2))) * 180/pi)
plot(time, wrapToPi(rpy_matlab(:,2)) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(3) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(3)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(rpy_(:,3)) * 180/pi)
plot(time, wrapToPi(data(:,ind_rpy(3))) * 180/pi)
plot(time, wrapToPi(rpy_matlab(:,3)) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

figure(5)
plot(time, bias_ * 180/pi), grid on
ylabel('Gyro Bias (deg/s)'), xlabel('Time (sec)')

figure(6)
subplot(211)
plot(time, data(:,ind_quat)), grid on
ylabel('Quaternion')
subplot(212)
plot(time, quat_), grid on
ylabel('Quaternion'), xlabel('Time (sec)')
