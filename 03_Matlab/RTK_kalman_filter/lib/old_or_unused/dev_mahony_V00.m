clc, clear variables
%%

data = readmatrix('putty_03.log');

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
kp = 1.0;
ki = kp * (2*pi* 0.8);
N = size(data, 1);
quat = eye(4,1);
bias = zeros(3,1);
p    = zeros(4,1);
Eva = [0 0 1].';
quat_ = zeros(N, 4);
bias_ = zeros(N, 3);

% https://ahrs.readthedocs.io/en/latest/filters/mahony.html
% dq/dt = 1/2 * q * p( Omega_y - b + kp * wmes) ~= ( q(k+1) - q(k) ) / Ts
% <-> q(k+1) = q(k) + Ts * 1/2 * q * p( Omega_y - b + kp * wmes)
% for i = 1:N
%     CBE = quat2CBE(quat);
% 
%     Bva  = CBE * Eva;
%     Bva_ = data(i, ind_acc).';
%     Bva_ = Bva_ ./ norm( Bva_ );
% 
%     w_mes = -ki * cross(Bva, Bva_);
%     
%     bias = bias + Ts * -ki * w_mes;
% 
%     p(2:4) = data(i,ind_gyro).' - bias + kp * w_mes;
% 
%     quat = quat + Ts * 1/2 * quat2QprodL( quat ) * p;
%     quat = quat ./ norm(quat);
% 
%     quat_(i,:) = quat.';
%     bias_(i,:) = bias.';
% end

% % betaflight (try)
% for i = 1:N
%     CEB = quat2CBE(quat).';
% 
%     Eva_ = CEB * data(i, ind_acc).';
%     Eva_ = Eva_ ./ norm( Eva_ );
% 
%     w_mes = -ki * cross(Eva, Eva_);
%     
%     bias = bias + Ts * -ki * w_mes;
% 
%     p(2:4) = CEB * data(i,ind_gyro).' - bias + kp * w_mes;
% 
%     quat = quat + Ts * 1/2 * quat2QprodL( quat ) * p;
%     quat = quat ./ norm(quat);
% 
%     quat_(i,:) = quat.';
%     bias_(i,:) = bias.';
% end

% betaflight
e = zeros(3,1);
integralFB = zeros(3,1);
for i = 1:N
    rMat = quat2CBE(quat).';

    a = data(i, ind_acc).';
    a = a ./ norm(a);
    ax = a(1); ay = a(2); az = a(3);
    e = [(ay * rMat(3,3) - az * rMat(3,2)); ...
         (az * rMat(3,1) - ax * rMat(3,3)); ...
         (ax * rMat(3,2) - ay * rMat(3,1))];

    integralFB = integralFB + ki * e * Ts;

    g = data(i, ind_gyro).' + kp * e + integralFB;
%     g = data(i, ind_gyro).';
%     g = g + kp * e + integralFB;
%     g = g * 1/2 * Ts;

    quat = quat + 1/2 * Ts * quat2QprodL( quat ) * [0; g];
    quat = quat ./ norm(quat);

%     buffer = quat;
%     quat = quat + [-buffer(2) * g(1) - buffer(3) * g(2) - buffer(4) * g(3); ...
%                    +buffer(1) * g(1) + buffer(3) * g(3) - buffer(4) * g(2); ...
%                    +buffer(1) * g(2) - buffer(2) * g(3) + buffer(4) * g(1); ...
%                    +buffer(1) * g(3) + buffer(2) * g(2) - buffer(3) * g(1)];
%     q.w += (-buffer.x * gx - buffer.y * gy - buffer.z * gz);
%     q.x += (+buffer.w * gx + buffer.y * gz - buffer.z * gy);
%     q.y += (+buffer.w * gy - buffer.x * gz + buffer.z * gx);
%     q.z += (+buffer.w * gz + buffer.x * gy - buffer.y * gx);

    quat_(i,:) = quat.';
    bias_(i,:) = integralFB.';
end

rpy_ = quat2rpy(quat_);

figure(44)
plot(time, cumtrapz(time, data(:,ind_gyro)) * 180/pi), grid on, hold on
plot(time, rpy_ * 180/pi), hold off
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')

figure(55)
plot(time, bias_ * 180/pi), grid on
ylabel('Gyro Bias (deg/s)'), xlabel('Time (sec)')
    
    
