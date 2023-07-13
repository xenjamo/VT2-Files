clc, clear all
%%
addpath('lib');
load mat_files/data_006.mat

A_mag = [ 0.9821680,  0.0000000,  0.0000000;
          0.0159774,  0.9866579,  0.0000000;
          0.0536142,  0.0066898,  1.0311741];
b_mag = [-0.4711387,  0.3239450, -0.3728509];
mag_calib = (mag - b_mag)*A_mag.';

tend = t(end-1);
t = t - t(1);
Ts = median(diff(t));

use_bf_data = 0;

%% get bf data
if use_bf_data
    mat_name = 'mat_files/data_20230609_bf_020';
    % mat_name = 'data_20230609_bf_021';
    % mat_name = 'data_20230609_bf_022';
    
    [ind, gyro_bf, acc_bf, mag_bf] = get_ind_and_synched_bf_data_20230609(gyro, mat_name);
    
    
    % uncomment if working with bf data
    t = t(ind,:);
    % gyro = gyro(ind,:);
    % accel = accel(ind,:);
    gyro = gyro_bf;
    accel = acc_bf;
    rtk_flags = rtk_flags(ind,:);
    relvelNED = relvelNED(ind,:);
    hpllh = hpllh(ind,:);
    K_pos = K_pos(ind,:);
    K_vel = K_vel(ind,:);
    mag_calib = mag_calib(ind,:);
    DOP = DOP(ind,:);
    numSV = numSV(ind,:);

end



%% orientation

% bessel
p = 2;         % pole at p rad/s
kp = 2 * p;
ki = kp^2 / 3;

para.kp = kp;
para.ki = ki;

rpy0 = 0 * [60, -60, 0] * 180/pi;
quat0 = rpy2quat(rpy0).';

if use_bf_data
    [quatRPY, biasRPY] = mahonyRPY(gyro_bf, acc_bf, mag_bf, para, Ts, quat0);
else
    [quatRPY, biasRPY] = mahonyRPY(gyro, accel, mag_calib, para, Ts, quat0);
    % [quatRP , biasRP ] = mahonyRP (gyro, accel, para, Ts, quat0);
end


% rpyRP  = quat2rpy(quatRP );
rpyRPY = quat2rpy(quatRPY);

% figure(6);
% subplot(211)
% plot(t, rpyRPY(:,1:2) * 180/pi), grid on
% title("rpyRPY")
% subplot(212)
% plot(t, rpyRPY(:,3)* 180/pi), grid on

% pmic: handle this in case you use the other gyro as source
rpyRPYint = cumtrapz(t, gyro) ;
figure(6);
subplot(211)
plot(t, rpyRPY(:,1:2) * 180/pi), grid on, hold on
plot(t, rpyRPYint(:,1:2)* 180/pi), hold off
title("rpyRPY")
subplot(212)
plot(t, unwrap(rpyRPY(:,3)) * 180/pi), grid on, hold on
plot(t, rpyRPYint(:,3)* 180/pi), hold off


% figure(61)
% subplot(411)
% plot(t, mag_calib); grid on;
% subplot(412)
% plot(t, mag_bf); grid on;
% subplot(413)
% plot(t, sqrt(sum(mag_calib.^2, 2))); grid on;
% subplot(414)
% plot(t, sqrt(sum(mag_bf.^2, 2))); grid on;

%% rotation
% acc_dev = 0.004879771823504,   0.004323229867404,   0.006564916157247
% 006 start angle 5.3965 rad / -2.70344
[m,n] = size(accel);
accel_M = zeros(m,n);
accel_ENU = zeros(m,n);
m_d = -3.36*pi/180; %magnetic declination
CEM = [cos(m_d) -sin(m_d) 0; sin(m_d) cos(m_d) 0; 0 0 1];

% vE = CEB * vB

for i = 1:m
    CMB = quat2rotm( quatRPY(i,:) ); % R
    accel_M(i,:) = CMB * accel(i,:)';
    accel_ENU(i,:) = CEM * accel_M(i,:)';
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

%% llh to NWU
lon = hpllh(:,1)*pi/180;
lat = hpllh(:,2)*pi/180;
h = hpllh(:,3);
% figure(7)
% plot3(lon , lat, h), grid on
% title('llh'), xlabel('lon'), ylabel('lat'), zlabel('h')


pos_ecef = transformWGS84ToECEF_R(lat, lon, h); % R stands for rad

% figure(8)
% plot3(pos_ecef(:,1) , pos_ecef(:,2) , pos_ecef(:,3)), grid on
% axis equal, title('pos ecef'), xlabel('x'), ylabel('y'), zlabel('z')

ind0 = 4;
pos_ecef_0 = transformWGS84ToECEF_R(lat(ind0), lon(ind0), h(ind0));
phi = lon(ind0);
la = lat(ind0);
R_ecefToLocal_0 = [ -sin(phi),          cos(phi),       0; ...
                    -cos(phi)*sin(la), -sin(la)*sin(phi), cos(la); ...
                     cos(la)*cos(phi),  cos(la)*sin(phi), sin(la)];
pos_nwu =  (pos_ecef - pos_ecef_0) * R_ecefToLocal_0.' * [0 -1 0;1 0 0; 0 0 1];

% figure(9)
% plot3(pos_nwu(:,1) , pos_nwu(:,2) , pos_nwu(:,3)), grid on,
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
figure(5)
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

Ac = [0 1 0; ...
     0 0 -1; ...
     0 0 0];
Bc = [0 1 0]';
Cc = [1 0 0; 0 1 0];
Dc = 0;

sys = ss(Ac,Bc,Cc,Dc);

Ad = eye(size(Ac)) + Ts * Ac;
Bd = Ts * Bc;
Cd = Cc;
Dd = Dc;

A = [Ad, zeros(3,6); zeros(3,3), Ad, zeros(3,3); zeros(3,6), Ad];
B = [Bd, zeros(3,2); zeros(3,1), Bd, zeros(3,1); zeros(3,2), Bd];
C = [Cd, zeros(size(Cd)), zeros(size(Cd));...
     zeros(size(Cd)), Cd, zeros(size(Cd));...
     zeros(size(Cd)), zeros(size(Cd)), Cd];
D = Dd;

x0_hat = [pos_nwu(1,1),  relvelNED(1,1) , 0,...
          pos_nwu(1,2), -relvelNED(1,2) , 0,...
          pos_nwu(1,3), -relvelNED(1,3) ,0]';
x_hat = zeros(m,9);
x = x0_hat; % zeros(3,1);
u = accel_M;
y = [pos_nwu(:,1),relvelNED(:,1),pos_nwu(:,2),-relvelNED(:,2),pos_nwu(:,3),-relvelNED(:,3)];
% y = [pos_nwu(:,1),pos_nwu(:,2),pos_nwu(:,3)];
e = zeros(size(y));
e_ = zeros(size(e(1,:)));

warn = zeros(m,1);
fault = zeros(m,1);
rtk_less = 0;

pos_nwu_kf = zeros(size(pos_nwu));
vel_nwu_kf = zeros(size(relvelNED));

P_ = ones(9,9)* 0.1;
P = zeros([size(x_hat),9]);

% var_gps = [hAcc.^2, hAcc.^2, vAcc.^2];
var_gps = [K_pos(:,1),K_pos(:,4),K_pos(:,6)];
var_vel = [K_vel(:,1),K_vel(:,4),K_vel(:,6)];
var_acc = diag([0 1 0 0 1 0 0 1 0] * 0.006564916157247*5);

% R = zeros(4824,6,6);
% R(:,1,1) = var_gps(:,1); R(:,2,2) = var_vel(:,1);
% R(:,3,3) = var_gps(:,2); R(:,4,4) = var_vel(:,2);
% R(:,5,5) = var_gps(:,3); R(:,6,6) = var_vel(:,3);
% R = [K_pos(:,1), K_vel(:,1), K_pos(:,4), K_vel(:,4), K_pos(:,6), K_vel(:,6)] / Ts;

Q = B * B.' * var_acc * Ts;
Q(3,3) = 1e1 * Ts; Q(6,6) = 1e1 * Ts; Q(9,9) = 1e1 * Ts;
rho = diag([1 3e-2 1 3e-2 1 3e-2]*0.1);
% k = 1583;
k = 3000;
% k = 1;
R = [K_pos(k,1), 0, K_pos(k,2), 0, K_pos(k,3), 0;...
     0 K_vel(k,1), 0 , K_vel(k,2), 0, K_vel(k,3);...
     K_pos(k,2), 0, K_pos(k,4), 0, K_pos(k,5), 0;...
     0 K_vel(k,2), 0 , K_vel(k,4), 0, K_vel(k,5);...
     K_pos(k,3), 0, K_pos(k,5), 0, K_pos(k,6), 0;...
     0 K_vel(k,3), 0 , K_vel(k,5), 0, K_vel(k,6)] / Ts * rho;
R_ = zeros(m,6);

K = dlqr(A', C', Q, R).'; % statische Lösung Kalman-Filter
sys_est = ss(A - K*C, [B, K], C, 0, Ts);

% eig_d = eig(sys_est.a);
% eig_c = log(eig_d) / Ts;
% fn = abs(eig_c) / 2 / pi
% Dn = cos(pi-angle(eig_c))
damp( sys_est )


for i = 2:m
    x = A * x + B * (u(i,:)'- [0 0 g].');
    P_ = A*P_*A' + Q;
    R_(i,:) = R_(i-1,:);
    if data_age(i) == 0
        R = [K_pos(i,1), 0, K_pos(i,2), 0, K_pos(i,3), 0;...
             0 K_vel(i,1), 0 , K_vel(i,2), 0, K_vel(i,3);...
             K_pos(i,2), 0, K_pos(i,4), 0, K_pos(i,5), 0;...
             0 K_vel(i,2), 0 , K_vel(i,4), 0, K_vel(i,5);...
             K_pos(i,3), 0, K_pos(i,5), 0, K_pos(i,6), 0;...
             0 K_vel(i,3), 0 , K_vel(i,5), 0, K_vel(i,6)] / Ts * rho;
        % R = diag([var_gps(i,1), var_vel(i,1), var_gps(i,2), var_vel(i,2), var_gps(i,3), var_vel(i,3)]) / Ts;
        % R = diag([var_gps(i,1), var_gps(i,2), var_gps(i,3)]) / Ts;
        R_(i,:) = eig(R)';
        e_ = y(i,:)' - C * x;
        S = C*P_*C' + R;
        K = P_* C' / S;
        x = x + K * e_;
        P_ = (eye(size(A))-K*C)*P_;
    end
    P(i,:,:) = P_;
    x_hat(i,:) = x';
    pos_nwu_kf(i,:) = [x(1), x(4), x(7)];
    vel_nwu_kf(i,:) = [x(2), x(5), x(8)];
    e(i,:) = e_';

    

    % pmic: carefull here, since we are not certain if everything on the
    % microcontroller side works as expected
    if data_age(i)*Ts > 0.4
        warn(i) = 1;
    end
    if data_age(i)*Ts > 0.8
        fault(i) = 1;
    end
    
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

    if max( eig(R) ) > 0.03
        fault(i) = 1;
    end



end

%%

figure(13)
subplot(311)
plot(t,[y(:,1),x_hat(:,1)]); grid on; title('pos n');
subplot(312);
plot(t,[y(:,3),x_hat(:,4)]); grid on; title('pos w');
subplot(313)
plot(t,[y(:,5),x_hat(:,7)]); grid on; title('pos u');

figure(1312)
subplot(311)
plot(t,[y(:,2),x_hat(:,2)]); grid on; title('vel n');
subplot(312);
plot(t,[y(:,4),x_hat(:,5)]); grid on; title('vel w');
subplot(313)
plot(t,[y(:,6),x_hat(:,8)]); grid on; title('vel u');

figure(131)
plot(t, x_hat(:,[3 6 9])), grid on, title('acc bias'); legend
% subplot(311)
% plot(t,x_hat(:,3)); grid on; title('acc bias n');
% subplot(312)
% plot(t,x_hat(:,6)); grid on; title('acc bias w');
% subplot(313)
% plot(t,x_hat(:,9)); grid on; title('acc bias u');

figure(132)
subplot(311)
plot(t,P(:,1,1)); hold on;
plot(t,P(:,2,2)); hold off;
grid on; title('covpos/covvel n'); legend({'covpos', 'covvel'}, 'location', 'best')
subplot(312)
plot(t,P(:,4,4)); hold on;
plot(t,P(:,5,5)); hold off;
grid on; title('covpos/covvel w'); legend({'covpos', 'covvel'}, 'location', 'best')
subplot(313)
plot(t,P(:,7,7)); hold on;
plot(t,P(:,8,8)); hold off;
grid on; title('covpos/covvel u'); legend({'covpos', 'covvel'}, 'location', 'best')

figure(133)
subplot(211)
plot(t, sqrt(var_gps)), grid on
grid on; title('covpos Ublox'); %legend({'covpos', 'covvel'}, 'location', 'best')
subplot(212)
plot(t, sqrt(var_vel)), grid on
grid on; title('covvel Ublox'); %legend({'covpos', 'covvel'}, 'location', 'best')
% subplot(311)
% plot(t,var_gps(:,1)); hold on;
% plot(t,var_vel(:,1)); hold off;
% grid on; title('covpos/covvel n Ublox'); legend({'covpos', 'covvel'}, 'location', 'best')
% subplot(312)
% plot(t,var_gps(:,2)); hold on;
% plot(t,var_vel(:,2)); hold off;
% grid on; title('covpos/covvel w Ublox'); legend({'covpos', 'covvel'}, 'location', 'best')
% subplot(313)
% plot(t,var_gps(:,3)); hold on;
% plot(t,var_vel(:,3)); hold off;
% grid on; title('covpos/covvel u Ublox'); legend({'covpos', 'covvel'}, 'location', 'best')

figure(134)
subplot(221)
plot(t, hAcc), grid on, ylabel('hAcc')
subplot(222)
plot(t, vAcc), grid on, ylabel('vAcc')
subplot(223)
plot(t, sAcc), grid on, ylabel('sAcc')
subplot(224)
plot(t, headAcc), grid on, ylabel('headAcc')




figure(19)
subplot(211)
plot(t,[e(:,1), e(:,3), e(:,5)]); grid on; title('error pos');
subplot(212)
plot(t,[e(:,2), e(:,4), e(:,6)]); grid on; title('error vel');
%% some stupid test

figure(1235)
subplot(211)
plot(t, [R_(:,1) R_(:,3) R_(:,5)]); grid on;
title('eigenwerte covpos')
subplot(212)
plot(t, [R_(:,2) R_(:,4) R_(:,6)]); grid on;
title('eigenwerte covvel')



%% hpllh
figure(14)
plot3(-pos_nwu(:,2),pos_nwu(:,1),pos_nwu(:,3), '.'); hold on;
plot3(-pos_nwu_kf(:,2),pos_nwu_kf(:,1),pos_nwu_kf(:,3),'-'); hold off;
axis equal, title("nwu"); xlabel('w'); ylabel('n'); zlabel('u')
grid on;

%%

figure(3)
plot3( pos_nwu_kf(:, 1), pos_nwu_kf(:, 2), pos_nwu_kf(:, 3) ), hold on
xlabel('x-Axis (m)'), ylabel('y-Axis (m)'), zlabel('z-Axis (m)')
for i = 1:50:m
    pos_i = pos_nwu_kf(i,:);
    CMB = quat2rotm ( quatRPY(i,:) ); % R
    arrow_length = 1;
    R_i_scaled = CMB * arrow_length;
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,1), R_i_scaled(2,1), R_i_scaled(3,1), 'LineWidth', 2, 'AutoScale', 'off', 'color', [1 0 0])
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,2), R_i_scaled(2,2), R_i_scaled(3,2), 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0.5 0])
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,3), R_i_scaled(2,3), R_i_scaled(3,3), 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0 1])
end
hold off, grid on, axis equal%, zlim([0 0.6])
set(findall(gcf, 'type', 'line'), 'linewidth', 1.5)

%%

figure(11);
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

figure(7)
subplot(311)
plot(t, carrSoln); grid on; title('carrSoln'); ylim([-0.2 2.2]);
subplot(312)
plot(t, warn); grid on; title('warn'); ylim([-0.2 1.2]);
subplot(313)
plot(t, fault); grid on; title('fault'); ylim([-0.2 1.2]);


% figure(15)
% subplot(211)
% plot(t, wrapToPi(-headMot))
% hold on;
% plot(t,rpyRPY(:,3));
% hold off
% title("headMot")
% grid on;
% subplot(212)
% plot(t, gSpeed); grid on
% title("gSpeed")
%% 

