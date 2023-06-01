clc
%%
load data.mat

data = [gyro, accel, mag_calib, t, quat, rpy];
tend = t(end-1);
time = data(:,10);
time = time - time(1);
Ts = median(diff(time));

ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
%%
w_gps = 1*2*pi;
w_rs = w_gps;

s = tf('s');
G_tp = w_gps/(s+w_gps);
G_hp = 1-G_tp;

G_tp_d = c2d(G_tp,Ts,'tustin');
G_hp_d = c2d(G_hp,Ts,'tustin');

figure(2)
bode(G_tp, G_hp)
%%

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

%%
figure(1);
t_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
nexttile(1,[2 1])
plot(time, unwrap(-headMot));grid on;hold on;
plot(time, unwrap(rpyRPY(:,3)))
plot(time, out.simout.data)
hold off;
title("heading");
legend(["headMot" "mahony Yaw" "fused"], "Location", "best")

gSpeed_acc = sqrt(sqrt(K_vel(:,1).^2 + K_vel(:,4).^2));
headAcc = gSpeed_acc./gSpeed; % dangerous may result in diffby0
headAcc = min(max(headAcc,0),pi);

nexttile(3)
plot(time, headAcc); grid on;




%%

b_a = 0;
A = [0 1 0; 0 0 -1; 0 0 b_a];
B = [0 1 0]';
C = [1 0 0];
D = 0;
