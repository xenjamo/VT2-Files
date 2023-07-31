clc, clear variables
%%

% data = readmatrix('putty_03.log');
% data = readmatrix('putty_04.log');
% data = readmatrix('putty_05.log');
data = readmatrix('putty_06.log');

time = data(:,end) * 1e-3;
time = time - time(1);
data = data(:,1:end-1);
Ts = median(diff(time));
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;

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

% Kp * (1 + 1/(Tn*s)) = kp + ki <->
% kp = Kp
% ki = Kp / Tn = kp / (1 / (2*pi*wn)) = kp * 2*pi*wn
kp = 10.0;
ki = kp * (2*pi* 1) * 0.5;
N = size(data, 1);
quat = eye(4,1);
bias = zeros(3,1);
p    = zeros(4,1);
Eva = [0 0 1].';
quat_ = zeros(N, 4);
bias_ = zeros(N, 3);

% https://ahrs.readthedocs.io/en/latest/filters/mahony.html
% adjusted so that it matches bf
% dq/dt = 1/2 * q * p( Omega_y - b + kp * wmes) ~= ( q(k+1) - q(k) ) / Ts
% <-> q(k+1) = q(k) + Ts * 1/2 * q * p( Omega_y - b + kp * wmes)
for i = 1:N
    CBE = quat2CBE(quat);

    Bva  = CBE * Eva;
    Bva_ = data(i, ind_acc).';
    Bva_ = Bva_ ./ norm( Bva_ );

    e = cross(Bva_, Bva);
    
    bias = bias + ki * e * Ts;

    p(2:4) = data(i,ind_gyro).' + bias + kp * e;

    quat = quat + Ts * 1/2 * quat2QprodL( quat ) * p;
    quat = quat ./ norm(quat);

    quat_(i,:) = quat.';
    bias_(i,:) = bias.';
end

% for i = 1:N
%     CEB = quat2CBE(quat).';
% 
%     Eva_ = CEB * data(i, ind_acc).';
%     Eva_ = Eva_ ./ norm( Eva_ );
% 
%     e = cross(Eva_, Eva);
%     
%     bias = bias + ki * e * Ts;
% 
%     p(2:4) = CEB * data(i,ind_gyro).' + bias + kp * e;
% 
%     quat = quat + Ts * 1/2 * quat2QprodL( quat ) * p;
%     quat = quat ./ norm(quat);
% 
%     quat_(i,:) = quat.';
%     bias_(i,:) = bias.';
% end

% % betaflight
% for i = 1:N
%     CBE = quat2CBE(quat);
%     
%     a = data(i, ind_acc).';
%     a = a ./ norm(a);
%     e = cross(a, CBE(:,3));
%     
%     bias = bias + ki * e * Ts;
% 
%     g = data(i, ind_gyro).' + kp * e + bias;
% 
%     quat = quat + 1/2 * Ts * quat2QprodL( quat ) * [0; g];
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
plot(time, wrapToPi(rpy_(:,1)) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(2) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(2)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(rpy_(:,2)) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
ax(3) = nexttile;
plot(time, wrapToPi(cumtrapz(time, data(:,ind_gyro(3)))) * 180/pi), grid on, hold on
plot(time, wrapToPi(rpy_(:,3)) * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

figure(5)
plot(time, bias_ * 180/pi), grid on
ylabel('Gyro Bias (deg/s)'), xlabel('Time (sec)')
    
    
