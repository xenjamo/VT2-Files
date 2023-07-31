clc, clear variables
addpath ..\99_fcn_bib\
%%

% addpath ../99_Matlab
% addpath ../21_Measurements/20220714

% data = read_bin_data('log_40805.bin');
load log_40805.mat %  save log_40805 data
Teval = [10 80];

data.time = data.ti; data = rmfield( data , 'ti');
data.Lidar = data.Lidar(:,3);

ind_eval = data.time >= Teval(1) & data.time <= Teval(2);
ti = data.time(ind_eval,:);
ti = ti - ti(1);

% without interpolation
Ts = median(diff(ti))
Ts = 1/50;
time = ti;
N = length(time);

figure(99)
plot(ti(1:end-1), diff(ti)), grid on

Tmean = 5.0;
ind0 = data.time > 1 & data.time < Tmean;
mag0 = median(data.mag(ind0,1:3)).';
Ymag0 = atan2(data.mag(ind0,1), data.mag(ind0,2));
Ymag0 = median(Ymag0);

% without interpolation
gyro    = data.gyr(  ind_eval,:);
acc     = data.acc(  ind_eval,:);
mag     = data.mag(  ind_eval,:);
Lidar   = data.Lidar(ind_eval,:);
OF      = data.OF(   ind_eval,:);
RS      = data.RS(   ind_eval,:);

RPY     = [data.est_RPY(ind_eval,:), data.est_Vxyz(ind_eval,:), data.est_xyz(ind_eval,:)];

FT      = data.cntrl_FT(ind_eval,:);
baro    = data.Baro(ind_eval,:);
RS_yaw  = data.RS_yaw(ind_eval,:);
RS      = data.RS(ind_eval,:);

%%

% Siehe:  ..\20_Literatur\TU_Chemnitz_Transformation_Matrices.pdf, S. 55&56

% extract RS orientation (quaternion -> rpy)
quat_RS = RS(:,4:7);  
ind_first = find(quat_RS(:,1) ~= 0, 1); % before RS is activated, the data contains zeros
quat_RS0 = quat_RS(ind_first,:).';
quat_RS0 = quat_RS0/norm(quat_RS0); % normalize
quat_RS0(2:end) = -quat_RS0(2:end); % conjugate
% rotate quaternion according to initial rotation and convert to rpy
RPY_RS = zeros(N,3);
for i = ind_first:N
    quat_RS_i = quat_RS(i,:)./sqrt(sum(quat_RS(i,:).^2, 2)); % normalize
    quat_RS(i,:) = (quat2QprodR(quat_RS0) * quat_RS_i.').';  % quaternion multiplication
    RPY_RS(i,:)  = quat2rpy(quat_RS(i,:));                   % extract rpy
end
RPY_RS  = unwrap(RPY_RS);

L0 = [0.1795, 0.0026, -0.0159];
% transform to cg
RS = RS(:,1:3);
for i = 1:length(RS)
    % CEB := R <-> CBE := R^T
    RS(i,:) = RS(i,:) - L0 * quat2CBE(quat_RS(i,:)); % RS = RS - R * L0
end

% % extract RS orientation (quaternion -> rpy)
% quat_RS = RS(:,4:7);  
% ind_first = find(quat_RS(:,1) ~= 0, 1); % before RS is activated, the data contains zeros
% quat_RS0 = quat_RS(ind_first,:).';
% quat_RS0 = quat_RS0/norm(quat_RS0); % normalize
% R_RS0 = quat2CEB(quat_RS0);
% % rotate quaternion according to initial rotation and convert to rpy
% RPY_RS = zeros(N,3);
% for i = ind_first:N
%     quat_RS_i = quat_RS(i,:)./sqrt(sum(quat_RS(i,:).^2, 2)); % normalize
%     R_RS_i = quat2CEB(quat_RS_i)*R_RS0.';
%     quat_RS(i,:) = rotm2quat(R_RS_i);
%     RPY_RS(i,:)  = quat2rpy(quat_RS(i,:));                   % extract rpy
% end
% RPY_RS  = unwrap(RPY_RS);
% 
% L0 = [0.1795, 0.0026, -0.0159];
% % transform to cg
% RS = RS(:,1:3);
% for i = 1:length(RS)
%     % CEB := R <-> CBE := R^T
%     RS(i,:) = RS(i,:) - L0 * quat2CBE(quat_RS(i,:)); % RS = RS - R * L0
% end

figure(100)
plot(time, [RPY(:,[1 2 3]), RPY_RS]*180/pi), grid on

% compensate OF with gyro and LiDAR
% -Z*(avg_flowx + wy)
% -Z*(avg_flowy - wx)
scale = 1.0;
vOF(:,1) = (OF(:,1) + scale*gyro(:,2)).*Lidar(:,1);
vOF(:,2) = (OF(:,2) - scale*gyro(:,1)).*Lidar(:,1);
% plot([vOF]), ylim([-4 4])

% compensate LiDAR to height
pzL = Lidar(:,1).*cos(RPY(:,1)).*cos(RPY(:,2));

% uncompensated with orientation (assume we flight always pointing in RealSense initial orientation)
pOF = cumtrapz(vOF(:,[1 2]))*Ts;
pOF = pOF + RS(1,[1 2]);

% define filter
wg = 2*pi*10.0;
Gf = c2d(tf(wg^2, [1 2*0.7*wg wg^2]), Ts, 'tustin');

% uncompensated with orientation (assume we fligh always pointing in RealSense initial orientation)
dRS = diff(RS)/Ts; dRS = [zeros(1,size(dRS,2));dRS];
dRS = filtfilt(Gf.num{1}, Gf.den{1}, dRS);

ddRS = diff(RS)/Ts; ddRS = [zeros(1,size(ddRS,2));ddRS];
ddRS = filtfilt(Gf.num{1}, Gf.den{1}, ddRS);

dLidar = diff(Lidar(:,1))/Ts; dLidar = [zeros(1,size(Lidar(:,1),2));dLidar];
dLidar = filtfilt(Gf.num{1}, Gf.den{1}, dLidar);

ddLidar = diff(dLidar)/Ts; ddLidar = [zeros(1,size(Lidar(:,1),2));ddLidar];
ddLidar = filtfilt(Gf.num{1}, Gf.den{1}, ddLidar);

dpzL = diff(pzL)/Ts; dpzL = [zeros(1,size(Lidar(:,1),2));dpzL];
dpzL = filtfilt(Gf.num{1}, Gf.den{1}, dpzL);

% dpzB = diff(baro)/Ts; dpzB = [zeros(1,size(Lidar(:,1),2));dpzB];
% dpzB = filtfilt(Gf.num{1}, Gf.den{1}, dpzB);

% transform to body velocity
for i = 1:length(dRS)
    dRS(i,:) = dRS(i,:) * (calcRn([0;0;1], RPY_RS(i,3))* ...
                           calcRn([0;1;0], RPY_RS(i,2))* ...
                           calcRn([1;0;0], RPY_RS(i,1)));
    % = ( (Rz_psi * Ry_the * Rx_phi).' * dRS(i,:).' ).'
    % =    dRS(i,:)   * (Rz_psi * Ry_the * Rx_phi)
end

% acc = -kv*v - (w x v) <-> v = 1/kv * (-acc - (w x v))
% - (w x v) = coriolis term
% vy*wz - vz*wy
% vz*wx - vx*wz
% vx*wy - vy*wx

% air resistance coefficient
kvx = 0.2;
kvy = 0.2;
kvz = 1;    % too noisy
vacc = -acc(:,[1 2 3]) + [zeros(N,2), 9.8*ones(N,1)];
vacc(:,1) = vacc(:,1) + 0*(vOF(:,2).*gyro(:,3) - dpzL.*gyro(:,2));
vacc(:,2) = vacc(:,2) + 0*(dpzL.*gyro(:,1) - vOF(:,1).*gyro(:,3));
vacc(:,1) = vacc(:,1)/kvx; vacc(:,2) = vacc(:,2)/kvy; vacc(:,3) = vacc(:,3)/kvz;

% gyr_acc_mag_FT_baro_RS_RScam_RPY_RS = [gyro(:,1:3), acc(:,1:3), mag(:,1:3), FT, baro, RS, RScam, RPY_RS];
% gyr_acc_mag_rspos_rsyaw_lidar = [gyro(:,1:3), acc(:,1:3), mag(:,1:3), RS, RS_yaw];
gyr_acc_mag_rspos_rsyaw_lidar = [gyro(:,1:3), acc(:,1:3), mag(:,1:3), RS, RPY_RS(:,3), Lidar];

Ymag = atan2(gyr_acc_mag_rspos_rsyaw_lidar(:,7), gyr_acc_mag_rspos_rsyaw_lidar(:,8));
Ymag = unwrap(Ymag - Ymag0);

ind = 1:3;
RPYint = cumtrapz(time, gyr_acc_mag_rspos_rsyaw_lidar(:,ind) - 0*mean(gyr_acc_mag_rspos_rsyaw_lidar(time < Tmean,ind)));

%%

data = gyr_acc_mag_rspos_rsyaw_lidar;
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
clear gyr_acc_mag

%% dont worry up to here

% parameters converted from ekf_rpyvp, 13.07.2022
para.Ts = Ts;
para.wg = 2*pi*0.1;
para.g  = 9.81;
para.kvx = 0.2;
para.kvy = 0.2;
% para.wa = 2*pi*0.7;

% 6 / 8 states
% [n_gyro; n_v; n_b_g; n_b_a]
var_fx = [0.01*[1 1], 10*[1 1], 0.02*[1 1], 1*[1 1]];
% var_fx = [0.1*[1 1], 10*[1 1], 0.02*[1 1], 1*[1 1]];
% [n_acc]
var_gy = [5.0*[1 1]];
rho = 0.1;

x_k = zeros(6,1);
[X_k, K_k, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, xi_k] = ... 
    ekf_RP_6states(data, x_k, para, var_fx, var_gy, rho);
% x_k = zeros(8,1);
% [X_k, K_k, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, xi_k] = ... 
%     ekf_RP_8states(data, x_k, para, var_fx, var_gy, rho); % THIS IS THE ACTUAL C++ IMPLEMENTATION

% Mahony
p = 2;         % pole at p rad/s
kp = 2 * [p, p, p].';
ki = kp.^2 / 3;
para.kp = kp;
para.ki = ki;
rpy0 = 0 * [60, 60, 60] * pi/180;
quat0 = rpy2quat(rpy0).';
[quatRP , biasRP ] = mahonyRP (data(:,ind_gyro), data(:,ind_acc), para, Ts, quat0);
[quatRPY, biasRPY] = mahonyRPY(data(:,ind_gyro), data(:,ind_acc), data(:,ind_mag), para, Ts, quat0);
rpyRP  = quat2rpy(quatRP );
rpyRPY = quat2rpy(quatRPY);

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

show_int = 1;

figure(3)
ax(1) = subplot(311);
stairs(time, [rpyRP(:,1), X_k(:,1), show_int*RPYint(:,1)]*180/pi), grid on
ylabel('Roll (deg)'), xlabel('Time (sec)')
ax(2) = subplot(312);
stairs(time, [rpyRP(:,2), X_k(:,2), show_int*RPYint(:,2)]*180/pi), grid on
ylabel('Pitch (deg)'), xlabel('Time (sec)')
ax(3) = subplot(313);
stairs(time, [Ymag, show_int*RPYint(:,3)]*180/pi), grid on
ylabel('Yaw (deg)')
xlabel('Time (s)')
linkaxes(ax, 'x'), clear ax
xlim([0 time(end)])

figure(4)
ax(1) = subplot(211);
stairs(time, [dRS(:,[1 2]), X_k(:,[3 4])]), grid on
legend('vx RS','vy RS','vx EKF','vy EKF')
ax(1) = subplot(212);
stairs(time, [X_k(:,[5 6])]), grid on
legend('b_gx','b_gy')
ylabel('Gyro Bias (rad/sec)')
xlabel('Time (s)')
linkaxes(ax, 'x'), clear ax
xlim([0 time(end)])

figure(5)
subplot(211)
semilogy(time, P_svd_k), grid on
subplot(212)
semilogy(time, xi_k), grid on
