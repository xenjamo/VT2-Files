clc, clear all
% '001.bin' % under bridge;
% '002.bin' % under bridge;
% '003.bin' % under bridge;
% '004.bin' % broken;
% '005.bin' % mag calib file;
% '006.bin' % general flight;
% '007.bin' % line walk;
% '008.bin' % line walk;
% '011.bin' % Roll Pitch Yaw;
% '013.bin' % speed;
% '014.bin' % speeeed;
% '015.bin' % XYZ;
% '016.bin' % XYZ;
% '017.bin' % chirp;
% '018.bin' % acro;
% '020.bin' % acro;
% '021.bin' % acro;
% '022.bin' % acro;
%%
addpath('lib');
load mat_files/data_014.mat

A_mag = [ 0.9821680,  0.0000000,  0.0000000;
          0.0159774,  0.9866579,  0.0000000;
          0.0536142,  0.0066898,  1.0311741];
b_mag = [-0.4711387,  0.3239450, -0.3728509];
mag_calib = (mag - b_mag) * A_mag.';

tend = t(end-1);
t = t - t(1);
Ts = median(diff(t));

%Teval = [47.5 52];
Teval = [0 89];
Tsave = [0 89];

%% piece data
use_bf_data = false;

% get bf data
if use_bf_data
    mat_name = 'mat_files/data_20230609_bf_014';
    % mat_name = 'data_20230609_bf_021'; mat_name = 'data_20230609_bf_022';
    [ind, gyro_bf, acc_bf, mag_bf] = get_ind_and_synched_bf_data_20230609(gyro, mat_name);
    gyro = gyro_bf;
    accel = acc_bf;
    mag_calib = mag_bf;
    t = t(ind,:);
end

ind = Teval(1) <= t & t <= Teval(2);
s_ind = Tsave(1) <= t & t <= Tsave(2);

t = t(ind,:);
rtk_flags = rtk_flags(ind,:);
relvelNED = relvelNED(ind,:);
hpllh = hpllh(ind,:);
K_pos = K_pos(ind,:);
K_vel = K_vel(ind,:);
DOP = DOP(ind,:);
numSV = numSV(ind,:);
hAcc = hAcc(ind,:);
headAcc = headAcc(ind,:);
vAcc = vAcc(ind,:);
sAcc = sAcc(ind,:);
gyro = gyro(ind,:);
accel = accel(ind,:);
mag_calib = mag_calib(ind,:);
headMot = headMot(ind,:);

t = t-t(1);

%% orientation

% p = 5;         % pole at p rad/s
% kp = p;
% ki = 0;

% % double real pole
% p = 1;         % pole at p rad/s
% kp = 2 * p;
% ki = kp^2 / 4;

% bessel
p = 2;         % pole at p rad/s
kp = 2 * p;
ki = kp^2 / 3;

para.kp = kp;
para.ki = ki;

rpy0 = 0 * [60, -60, 0] * 180/pi;
quat0 = rpy2quat(rpy0).';

[quatRPY, biasRPY] = mahonyRPY(gyro, accel, mag_calib, para, Ts, quat0);

% rpyRP  = quat2rpy(quatRP );
rpyRPY = quat2rpy(quatRPY);

% figure(3001); subplot(211) plot(t, rpyRPY(:,1:2) * 180/pi), grid on
% title("rpyRPY") subplot(212) plot(t, rpyRPY(:,3)* 180/pi), grid on

% pmic: handle this in case you use the other gyro as source
rpyRPYint = cumtrapz(t, gyro) ;

figure(3002);
subplot(211)
plot(t, rpyRPY(:,1:2) * 180/pi), grid on, hold on
plot(t, rpyRPYint(:,1:2)* 180/pi), hold off
title("rpyRPY")
subplot(212)
plot(t, unwrap(rpyRPY(:,3)) * 180/pi), grid on, hold on
plot(t, rpyRPYint(:,3)* 180/pi), hold off

figure(3033)
plot(t, -wrapToPi(headMot)), grid on; hold on;
plot(t, rpyRPY(:,3)), hold off;
%saveForLatex1d(downsample([-wrapToPi(headMot), rpyRPY(:,3), headAcc]*180/pi,20),downsample(t,20),'heading.csv', 1.2)

figure(3003);
subplot(411)
plot(t, mag_calib); grid on;
% subplot(412)
% plot(t,mag_bf); grid on;
subplot(413)
plot(t, sqrt(sum(mag_calib.^2, 2))); grid on;
% subplot(414)
% plot(t, sqrt(sum(mag_bf.^2, 2))); grid on;

%% rotation
% acc_dev = 0.004879771823504,   0.004323229867404,   0.006564916157247 006
% start angle 5.3965 rad / -2.70344
[m,n] = size(accel);
accel_M = zeros(m,n);
accel_NWU = zeros(m,n);
m_d = -3.36*pi/180; %magnetic declination
CEM = [cos(m_d) -sin(m_d) 0; sin(m_d) cos(m_d) 0; 0 0 1];

% vE = CEB * vB

for i = 1:m
    CMB = quat2rotm( quatRPY(i,:) ); % R
    accel_M(i,:) = CMB * accel(i,:)';
    accel_NWU(i,:) = CEM * accel_M(i,:)';
end

figure(3004);
subplot(211)
plot(t,accel); grid on;
title("accel")
subplot(212)
plot(t,accel_NWU); grid on;
title("accel global")

% figure(3005); plot(t,sqrt(sum(accel.^2, 2))-sqrt(sum(accel_ENU.^2, 2)));
% grid on; title("quantisation error")

%% llh to NWU
lon = hpllh(:,1)*pi/180;
lat = hpllh(:,2)*pi/180;
h = hpllh(:,3);
%saveForLatex2d(downsample(hpllh, 20),'csv/hpllh.csv',1.2);
% figure(3006) plot3(lon , lat, h), grid on title('llh'), xlabel('lon'),
% ylabel('lat'), zlabel('h')


pos_ecef = transformWGS84ToECEF_R(lat, lon, h); % R stands for rad
%saveForLatex3d(downsample(pos_ecef, 20),'csv/pos_ecef.csv',1.2);
figure(3007); plot3(pos_ecef(:,1) , pos_ecef(:,2) , pos_ecef(:,3)), grid on
axis equal, title('pos ecef'), xlabel('x'), ylabel('y'), zlabel('z')

ind0 = 4;
pos_ecef_0 = transformWGS84ToECEF_R(lat(ind0), lon(ind0), h(ind0));
phi = lon(ind0);
la = lat(ind0);
R_ecefToLocal_0 = [ -sin(phi),          cos(phi),       0; ...
                    -cos(phi)*sin(la), -sin(la)*sin(phi), cos(la); ...
                     cos(la)*cos(phi),  cos(la)*sin(phi), sin(la)];
pos_nwu =  (pos_ecef - pos_ecef_0) * R_ecefToLocal_0.' * [0 -1 0;1 0 0; 0 0 1];

%saveForLatex3d(downsample((pos_ecef - pos_ecef_0) * R_ecefToLocal_0.', 20),'csv/pos_nwu.csv',1.2);
% figure(3008) plot3(pos_nwu(:,1) , pos_nwu(:,2) , pos_nwu(:,3)), grid on,
% axis equal, title('pos ecef t'), xlabel('n'), ylabel('w'), zlabel('u')
%% data age (needs to be implemented on HW)
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
figure(3009)
plot(t,data_age*Ts); grid on;
title('time since last update')
%% RTK status

gnssFixOk = bitget(rtk_flags, 1);
diffSoln = bitget(rtk_flags, 2);
rtk_float = bitget(rtk_flags, 4);
rtk_fix = bitget(rtk_flags, 5);
carrSoln = bitor(rtk_float, 2*rtk_fix);


%% 3D test
% acc_dev = 0.004879771823504,   0.004323229867404,   0.006564916157247
g = 9.81;

Ac = [zeros(3,3),   eye(3,3), zeros(3,3);...
      zeros(3,3), zeros(3,3),  -eye(3,3);...
      zeros(3,3), zeros(3,3), zeros(3,3)];
Bc = [zeros(3,3),   eye(3,3), zeros(3,3)]';
Cc = [  eye(3,3), zeros(3,3), zeros(3,3);
      zeros(3,3),   eye(3,3), zeros(3,3)];
Dc = 0;

sys = ss(Ac,Bc,Cc,Dc);

Ad = eye(size(Ac)) + Ts * Ac;
Bd = Ts * Bc;
Cd = Cc;
Dd = Dc;

A = Ad;
B = Bd;
C = Cd;
D = Dd;
%%
x0_hat = [pos_nwu(1,:), relvelNED(1,:), zeros(1,3)]';
x_hat = zeros(m,9);
x = x0_hat; % zeros(3,1);
u = accel_M;
y = [pos_nwu,relvelNED.*[1 -1 -1]];
% y = [pos_nwu(:,1),pos_nwu(:,2),pos_nwu(:,3)];
e = zeros(size(y));
e_ = zeros(size(e(1,:)));

warn = zeros(m,1);
fault = zeros(m,1);
rtk_less = 0;

pos_nwu_kf = zeros(size(pos_nwu));
vel_nwu_kf = zeros(size(relvelNED));

P_ = ones(9,9)* 0;
P = zeros([size(x_hat),9]);

% var_gps = [hAcc.^2, hAcc.^2, vAcc.^2];
var_acc = diag([0 0 0 1 1 1 0 0 0] * 0.006564916157247*5);

% R = zeros(4824,6,6); R(:,1,1) = var_gps(:,1); R(:,2,2) = var_vel(:,1);
% R(:,3,3) = var_gps(:,2); R(:,4,4) = var_vel(:,2); R(:,5,5) =
% var_gps(:,3); R(:,6,6) = var_vel(:,3); R = [K_pos(:,1), K_vel(:,1),
% K_pos(:,4), K_vel(:,4), K_pos(:,6), K_vel(:,6)] / Ts;

psi = 10;
%psi = 1;
Q = B * B.' * var_acc * Ts * 10;
Q(7:9,7:9) = eye(3,3) * psi * Ts;

rho = diag([1 1 1 1e-2 1e-2 1e-2])*1e-1;
%rho = eye(6,6);

i = 1;
cov_gps = [K_pos(i,1), K_pos(i,2), K_pos(i,3);...
           K_pos(i,2), K_pos(i,4), K_pos(i,5);...
           K_pos(i,3), K_pos(i,5), K_pos(i,6)];
cov_vel = [K_vel(i,1), K_vel(i,2), K_vel(i,3);...
           K_vel(i,2), K_vel(i,4), K_vel(i,5);...
           K_vel(i,3), K_vel(i,5), K_vel(i,6)];
R = blkdiag(cov_gps, cov_vel)  / Ts * rho;
R_ = zeros(m,6);

K = dlqr(A', C', Q, R).'; % statische Lösung Kalman-Filter
sys_est = ss(A - K*C, [B, K], C, 0, Ts);

% eig_d = eig(sys_est.a); eig_c = log(eig_d) / Ts; fn = abs(eig_c) / 2 / pi
% Dn = cos(pi-angle(eig_c))
damp( sys_est )
figure(3999)
pzplot(sys_est)
grid on; axis equal;

%%

% A
% C
% Q, R

% x(k+1) = f(x,u) * Ts + x(k)
% P = A * P * A.' + W

% e = y - g(x,u)
% S^-1
% K = P*C.' * S^-1
% x = x + K*e
% P = (I - K*C) * P

%%

Ts_ds = 0; % downsamples Ts
for i = 2:m
    
    % update orientation
    CMB = quat2rotm( quatRPY(i,:) ); % R
    % accel_M(i,:) = CMB * accel(i,:)';
    % u(i,:) = accel_M(i,:);
    % accel_NWU(i,:) = CEM * accel_M(i,:)';
    % u(i,:) = accel(i,:);
    % g_body = CMB' * [0 0 g]';

    % update linearisation of system
    % A
    % C
    % Q (R is calculated below)
    Ac(4:6,7:9) = -CMB;
    Bc(4:6,:) = CMB;

    A = eye(size(Ac)) + Ts * Ac;
    B = Ts * Bc;
    C = Cc;
    D = Dc;

    Q = B * B.' * var_acc * Ts;
    Q(7:9,7:9) = eye(3,3)* psi * Ts;

    % x(k+1) = f(x,u) * Ts + x(k)
    % x = [p, v, b]^T
    p = x(1:3);
    v = x(4:6);
    b = x(7:9);
    f = [v; CMB * (accel(i,:).' - b) - [0 0 g].'; zeros(3,1)];
    x = f*Ts + x;
    
    % P = A * P * A.' + W
    % % x = A * x + (B * (u(i,:) - g_body')');
    P_ = A*P_*A' + Q;

    Ts_ds = Ts_ds + Ts;
    if mod(i,5) == 0 %data_age(i) == 0
        cov_gps = [K_pos(i,1), K_pos(i,2), K_pos(i,3);...
                   K_pos(i,2), K_pos(i,4), K_pos(i,5);...
                   K_pos(i,3), K_pos(i,5), K_pos(i,6)];
        cov_vel = [K_vel(i,1), K_vel(i,2), K_vel(i,3);...
                   K_vel(i,2), K_vel(i,4), K_vel(i,5);...
                   K_vel(i,3), K_vel(i,5), K_vel(i,6)];
        R = blkdiag(cov_gps, cov_vel)  / Ts_ds * rho;

        % e = y - g(x,u)
        % S^-1
        % K = P*C.' * S^-1
        % x = x + K*e
        % P = (I - K*C) * P

        e_ = y(i,:)' - C * x;
        S = C*P_*C' + R;
        K = P_* C' / S;
        x = x + K * e_;
        P_ = (eye(size(A)) - K*C) * P_;
        Ts_ds = 0;
    end
    P(i,:,:) = P_;
    x_hat(i,:) = x';
    pos_nwu_kf(i,:) = [x(1:3)];
    vel_nwu_kf(i,:) = [x(4:6)];
    e(i,:) = e_';

    

    % pmic: carefull here, since we are not certain if everything on the
    % microcontroller side works as expected
    % if data_age(i)*Ts > 0.4
    %     warn(i) = 1;
    % end
    % if data_age(i)*Ts > 0.8
    %     fault(i) = 1;
    % end
    
    % pmic: also evaluate hAcc, vAcc, etc.

    if carrSoln(i) == 0
        rtk_less = rtk_less +1;
        warn(i) = 1;
        if rtk_less*Ts > 5
            fault(i) = 1;
        end
    else
        rtk_less = 0;
    end

    if max([hAcc(i)^2, vAcc(i)^2]) > 0.03
        fault(i) = 1;
    end



end

%%

figure(3010)
subplot(311)
plot(t,[y(:,1),x_hat(:,1)]); grid on; title('pos n');
subplot(312);
plot(t,[y(:,2),x_hat(:,2)]); grid on; title('pos w');
subplot(313)
plot(t,[y(:,3),x_hat(:,3)]); grid on; title('pos u');

%saveForLatex1d([y(:,1:3),x_hat(:,1:3)],t,'phat.csv',1.2);

figure(3011)
subplot(311)
plot(t,[y(:,4),x_hat(:,4)]); grid on; title('vel n');
subplot(312);
plot(t,[y(:,5),x_hat(:,5)]); grid on; title('vel w');
subplot(313)
plot(t,[y(:,6),x_hat(:,6)]); grid on; title('vel u');

%saveForLatex1d([y(:,4:6),x_hat(:,4:6)],t,'vhat.csv',1.2);

%saveForLatex3d(downsample(y(:,1:3),5),'p_speed.csv',1.2);
%saveForLatex3d(x_hat(:,1:3),'phat_speed.csv',1.2);

%saveForLatex1d([y(:,1:3),x_hat(:,1:3)],t,'phat.csv',1.2);
%saveForLatex1d([y(:,4:6),x_hat(:,4:6)],t,'vhat.csv',1.2);

% saveForLatex1d([y(s_ind,1:3),x_hat(s_ind,1:3)],t(s_ind),'phat_bad_RTK.csv',1.2);
% saveForLatex1d([y(s_ind,4:6),x_hat(s_ind,4:6)],t(s_ind),'vnhat_bad_RTK.csv',1.2);


figure(3012)
plot(t, x_hat(:,7:9)), grid on, title('acc bias'); legend
%saveForLatex1d(downsample(x_hat(:,7:9),5),downsample(t,5),'bias_a.csv',1.2);
% subplot(311) plot(t,x_hat(:,3)); grid on; title('acc bias n');
% subplot(312) plot(t,x_hat(:,6)); grid on; title('acc bias w');
% subplot(313) plot(t,x_hat(:,9)); grid on; title('acc bias u');

figure(3013)
subplot(211)
plot(t,P(:,1,1)); hold on;
plot(t,P(:,2,2));
plot(t,P(:,3,3));hold off;
grid on; title('covpos');
subplot(212)
plot(t,P(:,4,4)); hold on;
plot(t,P(:,5,5)); 
plot(t,P(:,6,6)); hold off;
grid on; title('covvel');
%saveForLatex1d(downsample([P(:,1,1),P(:,2,2),P(:,3,3)],2),downsample(t,2),'covpos.csv',1.2)


figure(3014)
subplot(211)
plot(t, sqrt([K_pos(:,1),K_pos(:,4),K_pos(:,6)])), grid on
grid on; title('covpos Ublox'); %legend({'covpos', 'covvel'}, 'location', 'best')
subplot(212)
plot(t, sqrt([K_vel(:,1),K_vel(:,4),K_vel(:,6)])), grid on
grid on; title('covvel Ublox'); %legend({'covpos', 'covvel'}, 'location', 'best')
% subplot(311) plot(t,var_gps(:,1)); hold on; plot(t,var_vel(:,1)); hold
% off; grid on; title('covpos/covvel n Ublox'); legend({'covpos',
% 'covvel'}, 'location', 'best') subplot(312) plot(t,var_gps(:,2)); hold
% on; plot(t,var_vel(:,2)); hold off; grid on; title('covpos/covvel w
% Ublox'); legend({'covpos', 'covvel'}, 'location', 'best') subplot(313)
% plot(t,var_gps(:,3)); hold on; plot(t,var_vel(:,3)); hold off; grid on;
% title('covpos/covvel u Ublox'); legend({'covpos', 'covvel'}, 'location',
% 'best')
%saveForLatex1d(downsample([K_pos(:,1),K_pos(:,4),K_pos(:,6)],2),downsample(t,2),'covpos_ublox.csv',1.2)

figure(3015)
subplot(221)
plot(t, hAcc), grid on, ylabel('hAcc')
subplot(222)
plot(t, vAcc), grid on, ylabel('vAcc')
subplot(223)
plot(t, sAcc), grid on, ylabel('sAcc')
subplot(224)
plot(t, headAcc), grid on, ylabel('headAcc')




figure(3016)
subplot(211)
plot(t,[e(:,1:3)]); grid on; title('error pos');
subplot(212)
plot(t,[e(:,4:6)]); grid on; title('error vel');

%% hpllh
figure(3017)
plot3(-pos_nwu(:,2),pos_nwu(:,1),pos_nwu(:,3), '.'); hold on;
plot3(-pos_nwu_kf(:,2),pos_nwu_kf(:,1),pos_nwu_kf(:,3),'-'); hold off;
axis equal, title("nwu"); xlabel('w'); ylabel('n'); zlabel('u')
grid on;

%%

figure(3018)
plot3( pos_nwu_kf(:, 1), pos_nwu_kf(:, 2), pos_nwu_kf(:, 3) ), hold on
xlabel('x/N-Axis (m)'), ylabel('y/W-Axis (m)'), zlabel('z/U-Axis (m)')
for i = 1:50:m
    pos_i = pos_nwu_kf(i,:);
    CMB = quat2rotm ( quatRPY(i,:) ); % R
    arrow_length = 3;
    R_i_scaled = CMB * arrow_length;
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,1), R_i_scaled(2,1), R_i_scaled(3,1), 'LineWidth', 2, 'AutoScale', 'off', 'color', [1 0 0])
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,2), R_i_scaled(2,2), R_i_scaled(3,2), 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0.5 0])
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,3), R_i_scaled(2,3), R_i_scaled(3,3), 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0 1])
end
hold off, grid on, axis equal%, zlim([0 0.6])
view(-37.5, 30)
set(findall(gcf, 'type', 'line'), 'linewidth', 1.5)


%%

figure(3019);
subplot(211)
plot(t,DOP);
title("DOP")
legend(["gDOP" "pDOP" "tDOP" "vDOP" "hDOP" "nDOP" "eDOP"], "location", "northwest")
grid on;
subplot(212)
plot(t,numSV);
title("number of visible satellites")
ylim([0 max(numSV)+5])
grid on;

figure(3020)
subplot(311)
plot(t, carrSoln); grid on; title('carrSoln'); ylim([-0.2 2.2]);
subplot(312)
plot(t, warn); grid on; title('warn'); ylim([-0.2 1.2]);
subplot(313)
plot(t, fault); grid on; title('fault'); ylim([-0.2 1.2]);


% figure(3021) subplot(211) plot(t, wrapToPi(-headMot)) hold on;
% plot(t,rpyRPY(:,3)); hold off title("headMot") grid on; subplot(212)
% plot(t, gSpeed); grid on title("gSpeed")
%% only for 002
%saveForLatex3d(downsample(pos_nwu, 10),'bridgewalk.csv',1.4);
%xx = [hAcc, DOP(:,1), numSV, carrSoln, fault];
%saveForLatex1d(downsample(xx, 10),downsample(t, 10),'gps_health.csv',1.2);

% xx = [hAcc,vAcc, carrSoln];
% saveForLatex1d(xx(s_ind,:),t(s_ind),'fehler_flug.csv',1.2);
