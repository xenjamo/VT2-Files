function [X_k, K_k, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, RPY_k, xi_k, Kmat_k] ...
           = ekf_RPYVP_15states_no_baro_vel_wrt_earth_delayed_RSyaw_acc_inp(Z, x_k, para, var_fxn, var_gxn, rho)

global use_acc_for_pos_integration ax_past ay_past az_past
use_acc_for_pos_integration = false;
ax_past = 0;
ay_past = 0;
az_past = 0;

% acc for lin. model
Z0 = zeros(size(Z, 2));
Z0(6) = 9.81;

[Q_k, R_k] = calc_Q_R(x_k, Z0, para, var_fxn, var_gxn, rho);
F_k = calc_F_k(x_k, Z0, para);
H_k = calc_H_k(x_k, Z0, para);
P_k = 1 * para.scale_P0 * eye(15);

Q_k0 = Q_k;
R_k0 = R_k;
P_k0 = P_k;
F_k0 = F_k;
H_k0 = H_k;

X_k  = zeros(size(Z,1),length(x_k));
X_k(1,:) = x_k;
P_svd_k  = zeros(size(Z,1),15);
xi_k  = zeros(size(Z,1),1);
Kmat_k = zeros(size(Z,1),length(x_k)*4);

% delayed position estimates
nd_rspos = para.nd_rspos;
pos_past = zeros(nd_rspos+1, 3);
% pos_past = [px(k)   , py(k)   , pz(k)   ; ...
%             px(k-1) , py(k-1) , pz(k-1) ; ...
%                 .        .        .
%             px(k-nd), py(k-nd), pz(k-nd)]

% delayed angle psi estimates
nd_rspsi = para.nd_rspsi;
psi_past = zeros(nd_rspsi+1, 1);

for i = 1:size(X_k,1)
    
    F_k = calc_F_k(x_k, Z(i,:), para);
    H_k = calc_H_k(x_k, Z(i,:), para);
    [Q_k, R_k] = calc_Q_R(x_k, Z(i,:), para, var_fxn, var_gxn, rho);
    
    x_k = fxd(x_k, Z(i,:), para);
    P_k = F_k * P_k * F_k.' + Q_k;
        
    % update delayed angle psi estimation
    psi_past = [x_k([3]).'; psi_past(1:nd_rspsi,:)];
    
    % update delayed position estimation
    pos_past = [x_k([7 8 9]).'; pos_past(1:nd_rspos,:)];
        
    e = Z(i,[13 10 11 12]).' - gy(x_k, psi_past, nd_rspsi, pos_past, nd_rspos, Z(i,:), para);
    
    S_inv_k = eye(4)/(H_k*P_k*H_k.' + R_k);
    K_k = (P_k*H_k.')*S_inv_k;
            
    x_k = x_k + K_k*e;
    P_k = (eye(15) - K_k*H_k) * P_k;
    
    X_k(i,:) = x_k.';
    P_svd_k(i,:) = svd(P_k).';
    xi_k(i) = e.'*S_inv_k*e;
    Kmat_k(i,:) = K_k(:).';
    
end

RPY_k = X_k(:,1:3);

return

function [F_k, A] = calc_F_k(x, u, para)

phi = x(1);
theta = x(2);
psi = x(3);
vx = x(4);
vy = x(5);
vz = x(6);
px = x(7);
py = x(8);
pz = x(9);
b_gx = x(10);
b_gy = x(11);
b_gz = x(12);
b_ax = x(13);
b_ay = x(14);
b_az = x(15);

gyro_x = u(1);
gyro_y = u(2);
gyro_z = u(3);

wg_xy  = para.wg_xy;
wg_z  = para.wg_z;
g   = para.g;
kvx = para.kvx;
kvy = para.kvy;
wa_z  = para.wa_z;
wm  = para.wm;
wa_xy = para.wa_xy;

acc_x  = u(4);
acc_y  = u(5);
acc_z  = u(6);

A = [[                                     (sin(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta) - (cos(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta),                                              -(b_gz*cos(phi) - gyro_z*cos(phi) + b_gy*sin(phi) - gyro_y*sin(phi))/cos(theta)^2,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,     -1, -(sin(phi)*sin(theta))/cos(theta), -(cos(phi)*sin(theta))/cos(theta),                    0,                                                  0,                                                  0]
[                                                                                     cos(phi)*(b_gz - gyro_z) + sin(phi)*(b_gy - gyro_y),                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,                         -cos(phi),                          sin(phi),                    0,                                                  0,                                                  0]
[                                                           (sin(phi)*(b_gz - gyro_z))/cos(theta) - (cos(phi)*(b_gy - gyro_y))/cos(theta),                      - (cos(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta)^2 - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)^2,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,              -sin(phi)/cos(theta),              -cos(phi)/cos(theta),                    0,                                                  0,                                                  0]
[   (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_y - b_ay) + (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_z - b_az), cos(phi)*cos(psi)*cos(theta)*(acc_z - b_az) - cos(psi)*sin(theta)*(acc_x - b_ax) + cos(psi)*cos(theta)*sin(phi)*(acc_y - b_ay), (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) - (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_y - b_ay) - cos(theta)*sin(psi)*(acc_x - b_ax), 0, 0, 0, 0, 0, 0,      0,                                 0,                                 0, -cos(psi)*cos(theta),   cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta), - sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta)]
[ - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_y - b_ay) - (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_z - b_az), cos(phi)*cos(theta)*sin(psi)*(acc_z - b_az) - sin(psi)*sin(theta)*(acc_x - b_ax) + cos(theta)*sin(phi)*sin(psi)*(acc_y - b_ay), (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) - (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_y - b_ay) + cos(psi)*cos(theta)*(acc_x - b_ax), 0, 0, 0, 0, 0, 0,      0,                                 0,                                 0, -cos(theta)*sin(psi), - cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta),   cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)]
[                                                                 cos(phi)*cos(theta)*(acc_y - b_ay) - cos(theta)*sin(phi)*(acc_z - b_az),                          - cos(theta)*(acc_x - b_ax) - cos(phi)*sin(theta)*(acc_z - b_az) - sin(phi)*sin(theta)*(acc_y - b_ay),                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,                                 0,                                 0,           sin(theta),                               -cos(theta)*sin(phi),                               -cos(phi)*cos(theta)]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 1, 0, 0, 0, 0, 0,      0,                                 0,                                 0,                    0,                                                  0,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 1, 0, 0, 0, 0,      0,                                 0,                                 0,                    0,                                                  0,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 1, 0, 0, 0,      0,                                 0,                                 0,                    0,                                                  0,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0, -wg_xy,                                 0,                                 0,                    0,                                                  0,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,                            -wg_xy,                                 0,                    0,                                                  0,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,                                 0,                             -wg_z,                    0,                                                  0,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,                                 0,                                 0,               -wa_xy,                                                  0,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,                                 0,                                 0,                    0,                                             -wa_xy,                                                  0]
[                                                                                                                                       0,                                                                                                                              0,                                                                                                                                                                          0, 0, 0, 0, 0, 0, 0,      0,                                 0,                                 0,                    0,                                                  0,                                              -wa_z]];

F_k = eye(15) + para.Ts*A;

return

function H_k = calc_H_k(x, u, para)

phi = x(1);
theta = x(2);
psi = x(3);
vx = x(4);
vy = x(5);
vz = x(6);
px = x(7);
py = x(8);
pz = x(9);
b_gx = x(10);
b_gy = x(11);
b_gz = x(12);
b_ax = x(13);
b_ay = x(14);
b_az = x(15);

gyro_x = u(1);
gyro_y = u(2);
gyro_z = u(3);

wg_xy  = para.wg_xy;
wg_z  = para.wg_z;
g   = para.g;
kvx = para.kvx;
kvy = para.kvy;
wa_z  = para.wa_z;
wm  = para.wm;
wa_xy = para.wa_xy;


H_k = [[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
[ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
[ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]];

return

function [Qk, Rk] = calc_Q_R(x, u, para, var_fx, var_gx, rho)

phi = x(1);
theta = x(2);
psi = x(3);
vx = x(4);
vy = x(5);
vz = x(6);
px = x(7);
py = x(8);
pz = x(9);
b_gx = x(10);
b_gy = x(11);
b_gz = x(12);
b_ax = x(13);
b_ay = x(14);
b_az = x(15);

gyro_x = u(1);
gyro_y = u(2);
gyro_z = u(3);

wg_xy  = para.wg_xy;
wg_z  = para.wg_z;
g   = para.g;
kvx = para.kvx;
kvy = para.kvy;
wa_z  = para.wa_z;
wm  = para.wm;
wa_xy = para.wa_xy;


var_fx_0 = var_fx(1);
var_fx_1 = var_fx(2);
var_fx_2 = var_fx(3);
var_fx_3 = var_fx(4);
var_fx_4 = var_fx(5);
var_fx_5 = var_fx(6);
var_fx_6 = var_fx(7);
var_fx_7 = var_fx(8);
var_fx_8 = var_fx(9);
var_fx_9 = var_fx(10);
var_fx_10 = var_fx(11);
var_fx_11 = var_fx(12);
var_fx_12 = var_fx(13);
var_fx_13 = var_fx(14);
var_fx_14 = var_fx(15);

var_gy_0 = var_gx(1);
var_gy_1 = var_gx(2);
var_gy_2 = var_gx(3);
var_gy_3 = var_gx(4);

% n_acc_z: additive                                                                                                      0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0, var_fx_14]];

% n_acc_z: whitin CEB
Q = [[ (var_fx_0*cos(theta)^2 + var_fx_2*cos(phi)^2*sin(theta)^2 + var_fx_1*sin(phi)^2*sin(theta)^2)/cos(theta)^2, (cos(phi)*sin(phi)*sin(theta)*(var_fx_1 - var_fx_2))/cos(theta), -(sin(theta)*(var_fx_2 + var_fx_1*sin(phi)^2 - var_fx_2*sin(phi)^2))/(sin(theta)^2 - 1),                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0,         0,         0,         0,         0,         0]
[                                            (cos(phi)*sin(phi)*sin(theta)*(var_fx_1 - var_fx_2))/cos(theta),            var_fx_1 - var_fx_1*sin(phi)^2 + var_fx_2*sin(phi)^2,                                       (sin(2*phi)*(var_fx_1 - var_fx_2))/(2*cos(theta)),                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0,         0,         0,         0,         0,         0]
[                    -(sin(theta)*(var_fx_2 + var_fx_1*sin(phi)^2 - var_fx_2*sin(phi)^2))/(sin(theta)^2 - 1),               (sin(2*phi)*(var_fx_1 - var_fx_2))/(2*cos(theta)),              -(var_fx_2 + var_fx_1*sin(phi)^2 - var_fx_2*sin(phi)^2)/(sin(theta)^2 - 1),                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                          var_fx_4*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))^2 + var_fx_5*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))^2 + var_fx_3*cos(psi)^2*cos(theta)^2, var_fx_3*cos(psi)*cos(theta)^2*sin(psi) - var_fx_5*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - var_fx_4*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)), var_fx_5*cos(phi)*cos(theta)*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) - var_fx_3*cos(psi)*cos(theta)*sin(theta) - var_fx_4*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)),        0,        0,        0,        0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0, var_fx_3*cos(psi)*cos(theta)^2*sin(psi) - var_fx_5*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - var_fx_4*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)),                                                                                                          var_fx_4*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))^2 + var_fx_5*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))^2 + var_fx_3*cos(theta)^2*sin(psi)^2, var_fx_4*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) - var_fx_3*cos(theta)*sin(psi)*sin(theta) - var_fx_5*cos(phi)*cos(theta)*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)),        0,        0,        0,        0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                               var_fx_5*cos(phi)*cos(theta)*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) - var_fx_3*cos(psi)*cos(theta)*sin(theta) - var_fx_4*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)),                                                               var_fx_4*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) - var_fx_3*cos(theta)*sin(psi)*sin(theta) - var_fx_5*cos(phi)*cos(theta)*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)),                                                                                                                 var_fx_3*sin(theta)^2 + var_fx_5*cos(phi)^2*cos(theta)^2 + var_fx_4*cos(theta)^2*sin(phi)^2,        0,        0,        0,        0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0, var_fx_6,        0,        0,        0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0, var_fx_7,        0,        0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0, var_fx_8,        0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0, var_fx_9,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0, var_fx_10,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0,         0, var_fx_11,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0,         0,         0, var_fx_12,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0,         0,         0,         0, var_fx_13,         0]
[                                                                                                          0,                                                               0,                                                                                       0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                                                                                         0,                                                                                                                                                                                                           0,        0,        0,        0,        0,         0,         0,         0,         0, var_fx_14]];

% [~, A] = calc_F_k(x, u, para);

Qk = para.Ts * Q;
% Qk = para.Ts * ( Q + 1/2*( Q*A.' + A*Q )*para.Ts + 1/3*( A*Q*A.' )*para.Ts^2 );
% Qk = para.Ts * ( Q + 1/2*( Q*A.' + A*Q )*para.Ts + 1/3*( 1/2*(Q*A.'*A.' + A*A*Q) + A*Q*A.' )*para.Ts^2 + 1/8*(A*Q*A.'*A.' + A*A*Q*A.')*para.Ts^3 + 1/20*A*A*Q*A.'*A.'*para.Ts^4);
 
Rk = 1/para.Ts*[[ rho*var_gy_0,            0,            0,            0]
[            0, rho*var_gy_1,            0,            0]
[            0,            0, rho*var_gy_2,            0]
[            0,            0,            0, rho*var_gy_3]];

return

function x_k = fxd(x, u, para)

global use_acc_for_pos_integration ax_past ay_past az_past

phi = x(1);
theta = x(2);
psi = x(3);
vx = x(4);
vy = x(5);
vz = x(6);
px = x(7);
py = x(8);
pz = x(9);
b_gx = x(10);
b_gy = x(11);
b_gz = x(12);
b_ax = x(13);
b_ay = x(14);
b_az = x(15);

gyro_x = u(1);
gyro_y = u(2);
gyro_z = u(3);

wg_xy  = para.wg_xy;
wg_z  = para.wg_z;
g   = para.g;
kvx = para.kvx;
kvy = para.kvy;
wa_z  = para.wa_z;
wm  = para.wm;
wa_xy = para.wa_xy;

acc_x  = u(4);
acc_y  = u(5);
acc_z  = u(6);

if use_acc_for_pos_integration
%     x_k = x + para.Ts*[                                                        gyro_x - b_gx - (cos(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta) - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)
%                                                                                                                         sin(phi)*(b_gz - gyro_z) - cos(phi)*(b_gy - gyro_y)
%                                                                                             - (sin(phi)*(b_gy - gyro_y))/cos(theta) - (cos(phi)*(b_gz - gyro_z))/cos(theta)
%  (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) - (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_y - b_ay) + cos(psi)*cos(theta)*(acc_x - b_ax)
%  (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_y - b_ay) - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) + cos(theta)*sin(psi)*(acc_x - b_ax)
%                                                                     cos(phi)*cos(theta)*(acc_z - b_az) - sin(theta)*(acc_x - b_ax) - g + cos(theta)*sin(phi)*(acc_y - b_ay)
%                                                                                                                                                      ax_past*para.Ts/2 + vx % p(k) = p(k-1) + Ts*v(k-1) + Ts^2/2*a(k-1)
%                                                                                                                                                      ay_past*para.Ts/2 + vy
%                                                                                                                                                      az_past*para.Ts/2 + vz
%                                                                                                                                                                 -b_gx*wg_xy
%                                                                                                                                                                 -b_gy*wg_xy
%                                                                                                                                                                  -b_gz*wg_z
%                                                                                                                                                                 -b_ax*wa_xy
%                                                                                                                                                                 -b_ay*wa_xy
%                                                                                                                                                                  -b_az*wa_z];
% ax_past = (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) - (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_y - b_ay) + cos(psi)*cos(theta)*(acc_x - b_ax);
% ay_past = (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_y - b_ay) - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) + cos(theta)*sin(psi)*(acc_x - b_ax);
% az_past =                                                                    cos(phi)*cos(theta)*(acc_z - b_az) - sin(theta)*(acc_x - b_ax) - g + cos(theta)*sin(phi)*(acc_y - b_ay);
    a_k = [(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) - (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_y - b_ay) + cos(psi)*cos(theta)*(acc_x - b_ax)
           (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_y - b_ay) - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) + cos(theta)*sin(psi)*(acc_x - b_ax)                                                                          
                                                                             cos(phi)*cos(theta)*(acc_z - b_az) - sin(theta)*(acc_x - b_ax) - g + cos(theta)*sin(phi)*(acc_y - b_ay)];
    v_k = x(4:6) + para.Ts*a_k;
    p_k = x(7:9) + para.Ts*v_k + 0.5*para.Ts^2*a_k;
                                                                    
    x_k = x + para.Ts*[                                                        gyro_x - b_gx - (cos(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta) - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)
                                                                                                                            sin(phi)*(b_gz - gyro_z) - cos(phi)*(b_gy - gyro_y)
                                                                                                - (sin(phi)*(b_gy - gyro_y))/cos(theta) - (cos(phi)*(b_gz - gyro_z))/cos(theta)
     (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) - (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_y - b_ay) + cos(psi)*cos(theta)*(acc_x - b_ax)
     (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_y - b_ay) - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) + cos(theta)*sin(psi)*(acc_x - b_ax)
                                                                        cos(phi)*cos(theta)*(acc_z - b_az) - sin(theta)*(acc_x - b_ax) - g + cos(theta)*sin(phi)*(acc_y - b_ay)
                                                                                                                                                                             vx % p(k) = p(k-1) + Ts*v(k-1)
                                                                                                                                                                             vy
                                                                                                                                                                             vz
                                                                                                                                                                    -b_gx*wg_xy
                                                                                                                                                                    -b_gy*wg_xy
                                                                                                                                                                     -b_gz*wg_z
                                                                                                                                                                    -b_ax*wa_xy
                                                                                                                                                                    -b_ay*wa_xy
                                                                                                                                                                     -b_az*wa_z];
    x_k(7:9) = p_k;
else
    x_k = x + para.Ts*[                                                        gyro_x - b_gx - (cos(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta) - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)
                                                                                                                            sin(phi)*(b_gz - gyro_z) - cos(phi)*(b_gy - gyro_y)
                                                                                                - (sin(phi)*(b_gy - gyro_y))/cos(theta) - (cos(phi)*(b_gz - gyro_z))/cos(theta)
     (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) - (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_y - b_ay) + cos(psi)*cos(theta)*(acc_x - b_ax)
     (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_y - b_ay) - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) + cos(theta)*sin(psi)*(acc_x - b_ax)
                                                                        cos(phi)*cos(theta)*(acc_z - b_az) - sin(theta)*(acc_x - b_ax) - g + cos(theta)*sin(phi)*(acc_y - b_ay)
                                                                                                                                                                             vx % p(k) = p(k-1) + Ts*v(k-1)
                                                                                                                                                                             vy
                                                                                                                                                                             vz
                                                                                                                                                                    -b_gx*wg_xy
                                                                                                                                                                    -b_gy*wg_xy
                                                                                                                                                                     -b_gz*wg_z
                                                                                                                                                                    -b_ax*wa_xy
                                                                                                                                                                    -b_ay*wa_xy
                                                                                                                                                                     -b_az*wa_z];
end
return

function y_k = gy(x, psi_past, nd_rspsi, pos_past, nd_rspos, u, para)

phi = x(1);
theta = x(2);
psi = x(3);
vx = x(4);
vy = x(5);
vz = x(6);
px = x(7);
py = x(8);
pz = x(9);
b_gx = x(10);
b_gy = x(11);
b_gz = x(12);
b_ax = x(13);
b_ay = x(14);
b_az = x(15);

gyro_x = u(1);
gyro_y = u(2);
gyro_z = u(3);

wg_xy  = para.wg_xy;
wg_z  = para.wg_z;
g   = para.g;
kvx = para.kvx;
kvy = para.kvy;
wa_z  = para.wa_z;
wm  = para.wm;
wa_xy = para.wa_xy;

acc_x  = u(4);
acc_y  = u(5);
acc_z  = u(6);

% y_k = [psi
%   px
%   py
%   pz];

y_k = [psi_past(nd_rspsi+1, 1)
       pos_past(nd_rspos+1, 1)
       pos_past(nd_rspos+1, 2)
       pos_past(nd_rspos+1, 3)];

return