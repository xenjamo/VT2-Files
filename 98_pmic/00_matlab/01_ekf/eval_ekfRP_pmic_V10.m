clc, clear all
%%

% addpath ../99_Matlab
% addpath ../21_Measurements/20220714
% 
% data = read_bin_data('log_40805.bin');
load log_40805.mat
Teval = [8 85];

data.Lidar = data.Lidar(:,3);
data.RPY = data.est_RPY;
data.FT = data.cntrl_FT;

ind = data.ti >= Teval(1) & data.ti <= Teval(2);
ti = data.ti(ind,:);
ti = ti - ti(1);

Ts = 1/50;
time = (0:Ts:floor(ti(end)/Ts - 1)*Ts)';
N = length(time);
figure(99)
plot(ti(1:end-1), diff(ti)), grid on

Tmean = 1.0;
ind0 = data.ti > 1 & data.ti < Tmean;

Ymag0 = atan2(data.mag(ind0,1), data.mag(ind0,2));
Ymag0 = median(Ymag0);

% mag0(1) = mag0(1) - 1*0.03 + 0.2;
% mag0(2) = mag0(2) + 2*0.02 - 0.4;
% mag0(3) = mag0(3) - 3*0.01;

% data.mag(:,1) = data.mag(:,1) + 0.1;
% data.mag(:,2) = data.mag(:,2) + 0.1;

ipMethod = 'linear';
gyro  = interp1(ti, data.gyr(ind,:), time, ipMethod, 'extrap');
acc   = interp1(ti, data.acc(ind,:), time, ipMethod, 'extrap');
mag   = interp1(ti, data.mag(ind,:), time, ipMethod, 'extrap');
Lidar = interp1(ti, data.Lidar(ind,:), time, ipMethod, 'extrap');
OF    = interp1(ti, data.OF(ind,:), time, ipMethod, 'extrap');
RS    = interp1(ti, data.RS(ind,:), time, ipMethod, 'extrap');
RPY   = interp1(ti, data.RPY(ind,:), time, ipMethod, 'extrap');
FT    = interp1(ti, data.FT(ind,:), time, ipMethod, 'extrap');

% compensate OF with gyro and LiDAR
% -Z*(avg_flowx + wy)
% -Z*(avg_flowy - wx)
scale = 1.0;
vOF(:,1) = (OF(:,1) + scale*gyro(:,2)).*Lidar;
vOF(:,2) = (OF(:,2) - scale*gyro(:,1)).*Lidar;
% plot([vOF]), ylim([-4 4])

% compensate LiDAR to height
pz = Lidar.*cos(RPY(:,1)).*cos(RPY(:,2));

% uncompensated with orientation (assume we flight always pointing in RealSense initial orientation)
pOF = cumtrapz(vOF(:,[1 2]))*Ts;
pOF = pOF + RS(1,[1 2]);

% define filter
wg = 2*pi*2.0;
Gf = c2d(tf(wg^2, [1 2*0.7*wg wg^2]), Ts, 'tustin');

% uncompensated with orientation (assume we fligh always pointing in RealSense initial orientation)
dRS = diff(RS)/Ts; dRS = [zeros(1,size(dRS,2));dRS];
dRS = filtfilt(Gf.num{1}, Gf.den{1}, dRS);

ddRS = diff(RS)/Ts; ddRS = [zeros(1,size(ddRS,2));ddRS];
ddRS = filtfilt(Gf.num{1}, Gf.den{1}, ddRS);

dLidar = diff(Lidar)/Ts; dLidar = [zeros(1,size(Lidar,2));dLidar];
dLidar = filtfilt(Gf.num{1}, Gf.den{1}, dLidar);

ddLidar = diff(Lidar)/Ts; ddLidar = [zeros(1,size(Lidar,2));ddLidar];
ddLidar = filtfilt(Gf.num{1}, Gf.den{1}, ddLidar);

dpz = diff(pz)/Ts; dpz = [zeros(1,size(Lidar,2));dpz];
dpz = filtfilt(Gf.num{1}, Gf.den{1}, dpz);

% acc = -kv*v - (w x v) <-> v = 1/kv * (-acc - (w x v))
% - (w x v) = coriolis term
% vy*wz - vz*wy
% vz*wx - vx*wz
% vx*wy - vy*wx

% air resistance coefficient fg = 10, 5, 1, 0.5
kvx = 0.25; % robustfit(vacc(1/Ts:N-1/Ts,1), vOF(1/Ts:N-1/Ts,1))
kvy = 0.25; % robustfit(vacc(1/Ts:N-1/Ts,2), vOF(1/Ts:N-1/Ts,2))
kvz = 1;    % too noisy
vacc = -acc(:,[1 2 3]) + [zeros(N,2), 9.8*ones(N,1)];
vacc(:,1) = vacc(:,1) - 0*(vOF(:,2).*gyro(:,3) - dpz.*gyro(:,2));
vacc(:,2) = vacc(:,2) - 0*(dpz.*gyro(:,1) - vOF(:,1).*gyro(:,3));
vacc(:,1) = vacc(:,1)/kvx; vacc(:,2) = vacc(:,2)/kvy; vacc(:,3) = vacc(:,3)/kvz;
vacc = filtfilt(Gf.num{1}, Gf.den{1}, vacc);

gyr_acc_mag = [gyro(:,1:3), acc(:,1:3), mag(:,1:3)];

Ymag = atan2(gyr_acc_mag(:,7), gyr_acc_mag(:,8));
Ymag = unwrap(Ymag - Ymag0);

ind = 1:3;
RPYint = cumtrapz(time, gyr_acc_mag(:,ind) - 0*mean(gyr_acc_mag(time < Tmean,ind)));

% gyr_acc_mag = gyr_acc_mag(1:14,:);
% gyr_acc_mag(:,1) = [ 0.1,  0.5, -0.2,  0.2,  0.3,  0.2, -0.7, 0.5,  0.4, -0.3,  0.2,  0.8, -0.6,  0.2].';
% gyr_acc_mag(:,2) = [-0.1, -0.3,  0.1,  0.1,  0.1, -0.1, -0.9, 0.1,  0.5, -0.4, -0.2, -0.4,  0.5,  0.8].';
% gyr_acc_mag(:,3) = [ 0.5,  0.4, -0.3,  0.2,  0.8, -0.6,  0.2,-0.1, -0.3,  0.1,  0.1,  0.1, -0.1, -0.9].';
% gyr_acc_mag(:,4) = [ 0.3, -0.1,  0.7, -0.3, -0.9,  0.1,  0.2, 0.1,  0.5, -0.2,  0.2,  0.3,  0.2, -0.7].';
% gyr_acc_mag(:,5) = [ 0.1,  0.5, -0.4, -0.2, -0.4,  0.5,  0.8, 0.1, -0.2, -0.8,  0.2,  0.6,  0.4, -0.4].';
% gyr_acc_mag(:,7) = [ 0.1, -0.2, -0.8,  0.2,  0.6,  0.4, -0.4, 0.6, -0.4,  0.4,  0.7,  0.7,  0.5,  0.3].';
% gyr_acc_mag(:,8) = [ 0.6, -0.4,  0.4,  0.7,  0.7,  0.5,  0.3, 0.3, -0.1,  0.7, -0.3, -0.9,  0.1,  0.2].';

%%

% parameters converted from ekf_rpyvp, 13.07.2022
para.Ts = Ts;
para.wg = 2*pi*0.1*0;
para.g  = 9.81;
para.kvx = 0.2;
para.kvy = 0.2;
para.wa = 2*pi*0.7;

% 6 / 8 states
% [n_gyro; n_v; n_b_g; n_b_a]
var_fx = [0.01*[1 1], 10*[1 1], 0.02*[1 1], 1*[1 1]];
% var_fx = [1*[1 1], 1*[1 1], 1*[1 1], 1*[1 1]];
% [n_acc]
var_gy = [5.0*[1 1]];
rho = 0.1;

x_k = zeros(6,1);
[X_k, K_k, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, xi_k] = ... 
    ekf_RP_6states(gyr_acc_mag, x_k, para, var_fx, var_gy, rho); % THIS IS THE ACTUAL C++ IMPLEMENTATION
% x_k = zeros(8,1);
% [X_k, K_k, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, xi_k] = ... 
%     ekf_RP_8states(gyr_acc_mag, x_k, para, var_fx, var_gy, rho); % THIS IS THE ACTUAL C++ IMPLEMENTATION

% Statisches Kalman-Filter
format long
K_k_stat = dlqr(F_k0.', H_k0.', Q_k0, R_k0).'
format short

Ob = obsv(F_k0, H_k0);
rank(Ob)
cond(Ob)

format long
K_k
K_k_stat
norm(K_k - K_k_stat)
format short

eig_d = eig(F_k0 - K_k_stat*H_k0);
eig_c = log(eig_d)/Ts;
fn = abs(eig_c)/2/pi
Dn = cos(pi-angle(eig_c))

nx = size(F_k0,2);
sys_d = ss(F_k0 - K_k_stat*H_k0, [-Ts*[eye(2);zeros(nx-2,2)], -K_k_stat], eye(nx), zeros(nx,4), Ts);

figure(1)
sigma(sys_d), grid on, hold on
xlim([1e-5 1e2])
title('EKF Error Dynamics')

% figure(2)
% bodemag(sys_d([1 2 3 4],:)), grid on, hold on
% axis([1e-8 1e1 1e-6 1e2])
% title('EKF Error Dynamics')

%%

show_int = 1;

figure(3)
ax(1) = subplot(311);
stairs(time, [0*RPY(:,[1 2]), X_k(:,[1 2]), show_int*RPYint(:,[1 2])]*180/pi), grid on
ylabel('Angle (deg)')
ax(2) = subplot(312);
stairs(time, [Ymag, show_int*RPYint(:,3), RPY(:,3)]*180/pi), grid on
ylabel('Yaw-Angle (deg)')
ax(3) = subplot(313);
stairs(time, [dRS(:,[1 2]), X_k(:,[3 4]), vOF]), grid on
legend('vx RS','vy RS','vx EKF','vy EKF','vx OF','vy OF')
xlabel('Time (s)')
ylabel('Vel. (m/s)')
linkaxes(ax, 'x');
clear ax
xlim([0 time(end)])

figure(4)
ax(1) = subplot(311);
plot(time, [X_k(:,[3 4]), vOF]), grid on
ylim([-2 2])
ylabel('Vel. (m/s)')
legend('vx EKF','vy EKF','vx OF','vy OF')
ax(2) = subplot(312);
plot(time, X_k(:,[5 6])*1e3), grid on
ylabel('Gyro Bias (mrad/s)')
ax(3) = subplot(313);
semilogy(time, P_svd_k), grid on
xlabel('Time (s)')
ylabel('Sigma Values P')
linkaxes(ax, 'x');
clear ax
xlim([0 time(end)])

    figure(5)
    ax(1) = subplot(211);
    plot(time, X_k(:,[5 6])*1e3), grid on
    ylabel('Gyro Bias (mrad/s)')
    ax(2) = subplot(212);
if size(X_k,2) == 8
    plot(time, X_k(:,[7 8])*1e3), grid on
    xlabel('Time (s)')
    ylabel('Acc. Bias (mm/s^2)')
    linkaxes(ax, 'x');
    clear ax
    xlim([0 time(end)])
end

figure(6)
subplot(211)
semilogy(time, P_svd_k), grid on
subplot(212)
semilogy(time, xi_k), grid on

%%

% clc
% for i = 1:length(Q_k)
%     fprintf('%21.7ef, %14.7ef, %14.7ef, %14.7ef, %14.7ef, %14.7ef, %14.7ef, %14.7ef, \n', Q_k(i,:));
% end
% for i = 1:length(R_k)
%     fprintf('%21.7ef, %14.7ef, \n', R_k(i,:));
% end
