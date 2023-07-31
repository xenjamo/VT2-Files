clc, clear variables
addpath ..\99_fcn_bib\
%%

% addpath ../99_Matlab
% addpath ../21_Measurements/20220714

% data = read_bin_data('log_40805.bin');
load ..\log_40805.mat %  save log_40805 data
Teval = [0.02 inf];

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

% gyr_acc_mag_rspos_rsyaw_lidar = gyr_acc_mag_rspos_rsyaw_lidar(1:14,:);
% gyr_acc_mag_rspos_rsyaw_lidar(:,1)  = [  0.1,  0.5, -0.2,  0.2,  0.3,  0.2, -0.7,  0.5,  0.4, -0.3,  0.2,  0.8, -0.6,  0.2].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,2)  = [ -0.1, -0.3,  0.1,  0.1,  0.1, -0.1, -0.9,  0.1,  0.5, -0.4, -0.2, -0.4,  0.5,  0.8].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,3)  = [  0.5,  0.4, -0.3,  0.2,  0.8, -0.6,  0.2, -0.1, -0.3,  0.1,  0.1,  0.1, -0.1, -0.9].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,4)  = [  0.3, -0.1,  0.7, -0.3, -0.9,  0.1,  0.2,  0.1,  0.5, -0.2,  0.2,  0.3,  0.2, -0.7].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,5)  = [  0.1,  0.5, -0.4, -0.2, -0.4,  0.5,  0.8,  0.1, -0.2, -0.8,  0.2,  0.6,  0.4, -0.4].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,6)  = [ -0.1,  0.3,  0.7,  0.8,  0.9,  0.3, -0.3,  0.3,  0.5,  0.2,  0.1,  0.3, -0.7,  0.0].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,7)  = [  0.1, -0.2, -0.8,  0.2,  0.6,  0.4, -0.4,  0.6, -0.4,  0.4,  0.7,  0.7,  0.5,  0.3].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,8)  = [  0.6, -0.4,  0.4,  0.7,  0.7,  0.5,  0.3,  0.3, -0.1,  0.7, -0.3, -0.9,  0.1,  0.2].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,9)  = [  0.2,  0.1, -0.3,  0.3, -0.2,  0.4,  0.6, -0.2,  0.0, -0.1, -0.2, -0.1,  0.2,  0.3].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,10) = [  0.5, -0.7,  0.5,  0.6, -0.8, -0.2,  0.2,  0.5,  0.6,  0.6, -0.1,  0.4,  0.9,  0.7].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,11) = [  0.2, -0.0,  0.4, -0.2,  0.9, -0.5,  0.6,  0.6, -0.4,  0.7,  0.8, -0.2, -0.3, -0.2].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,12) = [ -0.7,  0.3,  0.2,  0.6, -0.7, -0.7,  0.5,  0.9,  0.9,  0.3,  0.3, -0.0, -0.5,  0.7].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,13) = [  0.3,  0.7,  0.9, -0.7, -3.2, -6.4, -3.2,  0.2, -0.1, -0.3,  0.0,  0.9,  0.7,  0.3].';
% gyr_acc_mag_rspos_rsyaw_lidar(:,13) = unwrap(gyr_acc_mag_rspos_rsyaw_lidar(:,13));

%%

% para.Ts = Ts;
% para.wg_xy = 2*pi*0.1*0;
% para.wg_z  = 2*pi*0.1*0;
% para.g  = 9.81;
% para.kvx = 0.150;
% para.kvy = 0.135;
% para.wm = 2*pi*0.8;
% para.wa_xy = 2*pi*0.7;
% para.wa_z = 2*pi*0;
% para.wrs_psi = 2*pi*0.8;
% para.wfil = 2*pi*15;
% para.nd_rspsi = 2;
% para.nd_rspos = 2;
% para.wrs = 2*pi*4.0;
% para.wrs_z = 2*pi*0.01;

para.Ts = Ts;
para.wg_xy = 2*pi*0.2
para.wg_z  = 2*pi*0.2;
para.g  = 9.81;
para.kvx = 0.2;
para.kvy = 0.2;
para.wm = 2*pi*0.8;
para.wa_xy = 2*pi*0.2;
para.wa_z = 2*pi*0;
para.wrs_psi = 2*pi*0.8;
para.wfil = 2*pi*15;
para.nd_rspsi = 2;
para.nd_rspos = 2;
para.wrs = 2*pi*4.0;
para.wrs_z = 2*pi*0.01;

para.scale_P0 = 1e2;

multpFigNr = 1;
var_fxn = [0.01*[1 1 1], 10*[1 1], 20, 1*[1 1], 1, 0.02*[1 1 10], 1*[1 1], 10, 4*[1 1], 1];
var_gyn = [5.0*[1 1], 1.0*[1 1], 3*[1 1 1]];
rho = 0.1;

x_k = zeros(17,1);
[X_k, Kk, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, RPY_k, xi_k, Kmat_k] = ... 
    ekf_RPYVP_17states_no_baro_vel_wrt_earth_delayed_01(gyr_acc_mag_rspos_rsyaw_lidar, x_k, para, var_fxn(1:17), var_gyn, mag0, rho);

% calculate OF with EKF altitude
vOF(:,1) = (OF(:,1) + scale*(gyro(:,2) - X_k(:,11))).*X_k(:,9)./(cos(X_k(:,1)).*cos(X_k(:,2)));
vOF(:,2) = (OF(:,2) - scale*(gyro(:,1) - X_k(:,10))).*X_k(:,9)./(cos(X_k(:,1)).*cos(X_k(:,2)));

% transform to body velocity
vel_X_k = X_k(:,[4 5 6]);
vel_RPY = RPY(:,[4 5 6]);
for i = 1:length(vel_X_k)
    vel_X_k(i,:) = vel_X_k(i,:) * (calcRn([0;0;1], RPY_k(i,3))* ...
                                   calcRn([0;1;0], RPY_k(i,2))* ...
                                   calcRn([1;0;0], RPY_k(i,1)));
    vel_RPY(i,:) = vel_RPY(i,:) * (calcRn([0;0;1], RPY(i,3))* ...
                                   calcRn([0;1;0], RPY(i,2))* ...
                                   calcRn([1;0;0], RPY(i,1)));                       
    % = ( (Rz_psi * Ry_the * Rx_phi).' * vel_k(i,:).' ).'
    % =    vek_k(i,:)   * (Rz_psi * Ry_the * Rx_phi)
end
X_k(:,[4 5 6]) = vel_X_k;
RPY(:,[4 5 6]) = vel_RPY;

% Statisches Kalman-Filter
K_k_stat = dlqr(F_k0.', H_k0.', Q_k0, R_k0).';

Ob = obsv(F_k0, H_k0);
rank(Ob)
cond(Ob)

% format long
% K_k
% K_k_stat
% norm(K_k - K_k_stat)
% format short

eig_d = eig(F_k0 - K_k_stat*H_k0);
eig_c = log(eig_d)/Ts;
fn = abs(eig_c)/2/pi
Dn = cos(pi-angle(eig_c))

if size(X_k, 2) == 17
    Bn0 = [eye(3), zeros(3,1); zeros(17-3,4)];
    Bn0(6,4) = 1;
    sys_d = ss(F_k0 - K_k_stat*H_k0, [-Bn0*Ts, -K_k_stat], eye(17), zeros(17,11), Ts);
elseif size(X_k, 2) == 15
    if size(H_k0, 1) ~= 4
        Bn0 = [eye(3), zeros(3,1); zeros(15-3,4)];
        Bn0(6,4) = 1;
    else
        Bn0 = zeros(15,6);
        Bn0(1:6,1:6) = eye(6);
    end
    sys_d = ss(F_k0 - K_k_stat*H_k0, [-Bn0*Ts, -K_k_stat], eye(15), zeros(15,10), Ts);
elseif size(X_k, 2) == 16
    Bn0 = [eye(3), zeros(3,1); zeros(16-3,4)];
    Bn0(6,4) = 1;
    sys_d = ss(F_k0 - K_k_stat*H_k0, [-Bn0*Ts, -K_k_stat], eye(16), zeros(16,10), Ts);
elseif size(X_k, 2) == 19
    Bn0 = zeros(19,6);
    Bn0(1:6,1:6) = eye(6);
    sys_d = ss(F_k0 - K_k_stat*H_k0, [-Bn0*Ts, -K_k_stat], eye(19), zeros(19,10), Ts);
end

figure(1)
sigma(sys_d), grid on, hold on
xlim([1e-4 1/2/Ts])

% return

%%

YLim_vxy = [-1 1]*2;

figure(creatFigNr(3, multpFigNr))
ax(1) = subplot(311);
% stairs(time, [RPY_RS(:,[1 2]), RPY_k(:,[1 2])]*180/pi), grid on
% legend('Roll RS', 'Pitch RS', 'Roll EKF', 'Pitch EKF', 'location', 'best')
stairs(time, [RPY_RS(:,[1 2]), RPY(:,[1 2]), RPY_k(:,[1 2])]*180/pi), grid on
ylabel('RP-Angle (deg)')
ax(2) = subplot(312);
% plot(time, [Ymag, RPYint(:,3), RPY_RS(:,3), RPY_k(:,3)]*180/pi), grid on
% legend('Magnetometer', 'Gyro-z integrated', 'RS', 'EKF', 'location', 'best')
plot(time, [Ymag, RPYint(:,3), RPY_RS(:,3), RPY(:,3), RPY_k(:,3)]*180/pi), grid on
ylabel('Y-Angle (deg)')
ax(3) = subplot(313);
% plot(time, [dRS(:,[1 2]), X_k(:,[4 5]), vOF]), grid on
% legend('Bvx RS', 'Bvy RS', 'Bvx EKF', 'Bvy EKF', 'Bvx OF', 'Bvy OF', 'location', 'best')
plot(time, [dRS(:,[1 2]), X_k(:,[4 5])]), grid on
legend('Bvx RS', 'Bvy RS', 'Bvx EKF', 'Bvy EKF', 'location', 'best')
xlabel('Time (s)')
ylabel('Vel. (m/s)')
ylim(YLim_vxy)
linkaxes(ax, 'x');
clear ax
xlim([0 time(end)])

figure(creatFigNr(4, multpFigNr))
ax(1) = subplot(311);
plot(time, [X_k(:,[4 5]), vOF]), grid on
ylabel('Vel. (m/s)')
legend('Bvx EKF', 'Bvy EKF', 'Bvx OF', 'Bvy OF', 'location', 'best')
ylim(YLim_vxy)
ax(2) = subplot(312);
plot(time, Kmat_k), grid on
ylim([-1 1]*0.5)
% plot(time, abs(Kmat_k)), grid on
% set(gca, 'YScale', 'log')
% ylim([1e-10 1e0])
ylabel('Kalman-Gain K_k')
ax(3) = subplot(313);
semilogy(time, P_svd_k), grid on
ylim([1e-3 1e4])
xlabel('Time (s)')
ylabel('Eigen Values P')
linkaxes(ax, 'x');
clear ax
xlim([0 time(end)])

figure(creatFigNr(5, multpFigNr))
ax(1) = subplot(311);
plot(time, X_k(:,10:12)*1e3), grid on
ylabel('Gyro Bias (mrad/s)')
legend('b_g_x', 'b_g_y', 'b_g_z', 'location', 'best')
ax(2) = subplot(312);
plot(time, X_k(:,13:15)), grid on
ylabel('Acc. Bias z (m/s^2)')
legend('b_a_x', 'b_a_y', 'b_a_z', 'location', 'best')
ax(3) = subplot(313);
if size(X_k, 2) == 17
    plot(time, X_k(:,16:17)), grid on
    ylabel('Mag. Bias (Gauss)')
elseif size(X_k, 2) == 15 || size(X_k, 2) == 19
    plot(time, zeros(size(time, 1), 2)), grid on
    ylabel('Mag. Bias (Gauss)')
elseif size(X_k, 2) == 16
    plot(time, X_k(:,16)), grid on
    ylabel('RS z Bias (deg)')
%     plot(time, X_k(:,16)*180/pi), grid on
%     ylabel('RS Yaw Bias (deg)')
end
xlabel('Time (s)')
linkaxes(ax, 'x');
clear ax
xlim([0 time(end)])

figure(creatFigNr(6, multpFigNr))
ax(1) = subplot(231);
plot(time, [dRS(:,1), RPY(:,4), X_k(:,4)]), grid on
ylim(YLim_vxy)
ylabel('Bvx (m/s)')
ax(2) = subplot(234);
plot(time, [RS(:,1), RPY(:,7), X_k(:,7)]), grid on
ylabel('Ipx (m)')
xlabel('Time (s)')
ax(3) = subplot(232);
plot(time, [dRS(:,2), RPY(:,5), X_k(:,5)]), grid on
ylim(YLim_vxy)
ylabel('Bvy (m/s)')
ax(4) = subplot(235);
% plot(time, [RS(:,2), X_k(:,8)]), grid on
plot(time, [RS(:,2), RPY(:,8), X_k(:,8)]), grid on
ylabel('Ipy (m)')
xlabel('Time (s)')
ax(5) = subplot(233);
% plot(time, [0*dpzB, dpzL, dRS(:,3), X_k(:,6)]), grid on
plot(time, [dRS(:,3), RPY(:,6), X_k(:,6)]), grid on
ylabel('Bvz (m/s)')
ylim([-1 1])
ax(6) = subplot(236);
% plot(time, [0*baro, pzL, RS(:,3), X_k(:,9)]), grid on
% plot(time, [RS(:,3), X_k(:,9)]), grid on
plot(time, [RS(:,3), RPY(:,9), X_k(:,9)]), grid on
ylabel('Ipz (m)')
xlabel('Time (s)')
legend('RS', 'EKF online', 'EKF', 'location', 'best')
linkaxes(ax, 'x');
clear ax
xlim([0 time(end)])

figure(creatFigNr(7, multpFigNr))
ax(1) = subplot(211);
semilogy(time, P_svd_k), grid on
ylim([1e-3 1e4])
ylabel('Eigen Values P')
ax(2) = subplot(212);
semilogy(time, xi_k), grid on
linkaxes(ax, 'x');
clear ax
xlim([0 time(end)])

figure(creatFigNr(8, multpFigNr))
plot(time, baro(:,1), '.'), grid on, hold on
plot(time, [pzL, RS(:,3), RPY(:,9), X_k(:,9)], '.', 'Linewidth', 2), hold off
ylabel('Ipz (m)')
xlabel('Time (s)')
legend('Baro', 'Lidar', 'RS', 'EKF online', 'EKF', 'location', 'best')
xlim([0 time(end)])
