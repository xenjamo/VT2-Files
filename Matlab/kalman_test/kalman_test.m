clc, clear all
%%
addpath('lib');
load mat_files/data_006.mat

A_mag = [ 0.9821680,  0.0000000,  0.0000000;
          0.0159774,  0.9866579,  0.0000000;
          0.0536142,  0.0066898,  1.0311741];
b_mag = [-0.4711387,  0.3239450, -0.3728509];
mag_calib = (mag - b_mag)*A_mag.';

data = [gyro, accel, mag_calib, t, quat, rpy];
tend = t(end-1);
time = data(:,10);
time = time - time(1);
Ts = median(diff(time));

ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
%% komplementär filter
w_gps = 1*2*pi;
w_rs = w_gps;

s = tf('s');
G_tp = w_gps/(s+w_gps);
G_hp = 1-G_tp;

G_tp_d = c2d(G_tp,Ts,'tustin');
G_hp_d = c2d(G_hp,Ts,'tustin');

% figure(2)
% bode(G_tp, G_hp)
%% orientation

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

figure(6);
plot(t,rpyRPY);
title("rpyRPY")
grid on;

%% rotation
% acc_dev = 0.004879771823504,   0.004323229867404,   0.006564916157247
% 006 start angle 5.3965 rad / -2.70344
[m,n] = size(accel);
accel_M = zeros(m,n);
accel_ENU = zeros(m,n);
m_d = -3.36*pi/180; %magnetic declination
CEM = [cos(m_d) -sin(m_d) 0; sin(m_d) cos(m_d) 0; 0 0 1];

for i = 1:m
    CMB = rotmat(quaternion(quatRP(i,:)),'frame');
    accel_M(i,:) = CMB*accel(i,:)';
    accel_ENU(i,:) = CEM*accel_M(i,:)';

end

figure(2);
subplot(211)
plot(t,accel); grid on;
title("accel")
subplot(212)
plot(t,accel_ENU); grid on;
title("accel global")

% figure(25);
% plot(t,sqrt(sum(accel.^2, 2))-sqrt(sum(accel_ENU.^2, 2))); grid on;
% title("quantisation error")
%% data age
dpos = [0 0 0;diff(relposNED)];
data_age = zeros(m,1);
age = 0;
for i = 1:m
    if dpos(i,:) == 0
        age = age +1;
    else
        age = 0;
    end
    data_age(i) = age;
end
figure(5)
plot(t,data_age*Ts); grid on;
title('time since last update')
%% 1d test
% acc_dev = 0.004879771823504,   0.004323229867404,   0.006564916157247
g = 9.81;

Ac = [0 1 0; ...
     0 0 -1; ...
     0 0 0];
Bc = [0 1 0]';
Cc = [1 0 0;0 1 0];
Dc = 0;

sys = ss(Ac,Bc,Cc,Dc);

Ad = eye(size(Ac)) + Ts * Ac;
Bd = Ts * Bc;
Cd = Cc;
Dd = Dc;

x0_hat = [hpllh(1,3), -relvelNED(1,3) ,0]';
x_hat = zeros(m,3);
x = x0_hat; % zeros(3,1);
u = accel_ENU(:,3);
y = [hpllh(:,3),-relvelNED(:,3)];

P_ = eye(3,3);
P = zeros([size(x_hat),3]);

var_gps = K_pos(:,6);
var_vel = K_vel(:,6);
var_acc = 0.006564916157247;

R = [var_gps, var_vel] / Ts;
Q = Bc * Bc.' * var_acc * Ts;
Q(3,3) = 1 * Ts;

for i = 1:m
    x = Ad * x + Bd * (u(i)-g);
    P_ = Ad*P_*Ad' + Q;
    
    if data_age(i) == 0
        e = y(i,:)' - Cd * x;
        S = Cd*P_*Cd' + R(i);
        K = P_* Cd' / S;
        x = x + K * e;
        P_ = (eye(size(Ad))-K*Cd)*P_;
    end
    P(i,:,:) = P_;
    x_hat(i,:) = x';
    
end
figure(13)
subplot(311)
plot(t,[y(:,1),x_hat(:,1)]); grid on; title('pos');
subplot(312);
plot(t,[y(:,2),x_hat(:,2)]); grid on; title('vel');
subplot(313)
plot(t,u); grid on; title('acc');

figure(131)
subplot(211)
plot(t,x_hat(:,3)); grid on; title('acc bias');
subplot(212)
plot(t,P(:,1,1)); hold on;
plot(t,P(:,2,2)); hold off;
grid on; title('covpos/covvel'); legend({'covpos', 'covvel'}, 'location', 'best')
%% hpposllh
figure(14)
plot3(hpllh(:,1),hpllh(:,2),hpllh(:,3), '.');
title("hpllh")
grid on;

figure(15)
subplot(211)
plot(t, wrapToPi(-headMot))
hold on;
plot(t,rpyRPY(:,3));
hold off
title("headMot")
grid on;
subplot(212)
plot(t, gSpeed); grid on
title("gSpeed")
