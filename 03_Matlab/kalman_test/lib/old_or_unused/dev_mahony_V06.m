clc, clear variables
%%

% data = readmatrix('putty_09.log');
% data = readmatrix('putty_10.log');
% data = readmatrix('putty_11.log');
% data = readmatrix('putty_12.log');
data = readmatrix('putty_13.log');

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
p = [2, 2, 2].';         % pole at p rad/s
kp = 2 * p;
ki = kp.^2 / 3 * 0;

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

rpy0 = 0 * [60, -60, 0] * 180/pi;
quat = rpy2quat(rpy0).';

% betaflight (rpy)
for i = 1:N

    % h = (hx; hy;  0) - measured mag field vector in EF (assuming Z-component is zero)
    % b = (bx;  0;  0) - reference mag field vector heading due North in EF (assuming Z-component is zero)
    mag_n = data(i, ind_mag).';
    mag_n = mag_n ./ norm(mag_n);

    CBE = quat2CBE(quat);
    CEB = CBE.';
    mu = [0 0 1].';
    Cmu = eye(3) - mu*mu.';
    h = Cmu.' * CEB * mag_n; % Cmu.' * CEB = (CEB.' * Cmu).' = = (CBE * Cmu).'

%     % by normalising here you adjust implicit the weigth
%     if i == 1
%         h_norm_init = norm(h);
%     end
%     b = [1; 0; 0];
%     e = CBE * cross(h * h_norm_init / norm(h), b);
    b = [1; 0; 0];
    e = CBE * cross(h/norm(h), b);
%     % by normalising here you adjust implicit the weigth
%     b = [norm(h); 0; 0];
%     e = CBE * cross(h, b);

    g_n = CEB(3,:).';
%     g_n = [         2*(quat(2)*quat(4) - quat(1)*quat(3)); ...
%                     2*(quat(3)*quat(4) + quat(1)*quat(2)); ...
%             quat(1)^2 - quat(2)^2 - quat(3)^2 + quat(4)^2];

    acc_n = data(i, ind_acc).';
    acc_n = acc_n ./ norm(acc_n);
    e = e + cross(acc_n, g_n);
        
    bias = bias + ki .* e * Ts;

    Q = [[-quat(2), -quat(3), -quat(4)]; ...
         [ quat(1), -quat(4),  quat(3)]; ...
         [ quat(4),  quat(1), -quat(2)]; ...
         [-quat(3),  quat(2),  quat(1)]];

    dquat = Ts * 0.5 * Q * ( data(i,ind_gyro).' + bias + kp .* e );
    quat = quat + dquat;
    quat = quat ./ norm(quat);

    quat_(i,:) = quat.';
    bias_(i,:) = bias.';
end

% % betaflight (rpy)
% for i = 1:N
% 
%     mag_n = data(i, ind_mag).';
%     mag_n = mag_n ./ norm(mag_n);
% 
%     rMat = quat2CBE(quat).';
%     hx = rMat(1,:) * mag_n;
%     hy = rMat(2,:) * mag_n;
%     bx = sqrt(hx * hx + hy * hy);
%     ez_ef = -(hy * bx);
%     e = rMat(3,:).' * ez_ef; % rMat(3,:) = CEB(3,:) = CBE(:,3)
%     
%     g_n = [         2*(quat(2)*quat(4) - quat(1)*quat(3)); ...
%                     2*(quat(3)*quat(4) + quat(1)*quat(2)); ...
%             quat(1)^2 - quat(2)^2 - quat(3)^2 + quat(4)^2];
% 
%     acc_n = data(i, ind_acc).';
%     acc_n = acc_n ./ norm(acc_n);
%     e = e + cross(acc_n, g_n);
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

rpy_ = quat2rpy(quat_);

figure(4)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(1)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(rpy_(:,1)) * 180/pi)
plot(time, wrapToPi(data(:,ind_rpy(1))) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(2) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(2)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(rpy_(:,2)) * 180/pi)
plot(time, wrapToPi(data(:,ind_rpy(2))) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(3) = nexttile;
plot(time, unwrap(cumtrapz(time, data(:,ind_gyro(3)))) * 180/pi), grid on, hold on
plot(time, unwrap(rpy_(:,3)) * 180/pi)
plot(time, unwrap(data(:,ind_rpy(3))) * 180/pi), hold off
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
