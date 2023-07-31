function [X_k, K_k, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, xi_k] ...
                   = ekf_RP_8states(Z, x_k, para, var_fx, var_gy, rho)
               
[Q_k, R_k] = calc_Q_R_ekf_test(zeros(8), zeros(2), para, var_fx, var_gy, rho);
F_k = calc_F_k(zeros(8), zeros(2), para);
H_k = calc_H_k(zeros(8), zeros(2), para);
P_k = 1 * 1e2 * eye(8);

Q_k0 = Q_k;
R_k0 = R_k;
P_k0 = P_k;
F_k0 = F_k;
H_k0 = H_k;

X_k  = zeros(size(Z,1), length(x_k));
X_k(1,:) = x_k;
P_svd_k  = zeros(size(Z,1), 8);
xi_k  = zeros(size(Z,1), 1);

for i = 1:size(X_k,1)
    
    F_k = calc_F_k(x_k, Z(i,:), para);
    H_k = calc_H_k(x_k, Z(i,:), para);
    [Q_k, R_k] = calc_Q_R_ekf_test(x_k, Z(i,:), para, var_fx, var_gy, rho);
    
    x_k = fxd(x_k, Z(i,:), para);
    P_k = F_k * P_k * F_k.' + Q_k;
    e = Z(i,[4 5]).' - gy(x_k, Z(i,:), para);
    
    S_inv_k = eye(2)/(H_k*P_k*H_k.' + R_k);
    K_k = (P_k*H_k.')*S_inv_k;
            
    x_k = x_k + K_k*e;
    P_k = (eye(8) - K_k*H_k) * P_k;
    
    X_k(i,:) = x_k.';
    P_svd_k(i,:) = svd(P_k).';
    xi_k(i) = e.'*S_inv_k*e;
    
end

return

function [F_k, A] = calc_F_k(x, u, para)

phi   = x(1);
theta = x(2);
vx    = x(3);
vy    = x(4);
b_gx  = x(5);
b_gy  = x(6);
b_ax  = x(7);
b_ay  = x(8);

gyro_x = u(1);
gyro_y = u(2);

wg  = para.wg;
g   = para.g;
kvx = para.kvx;
kvy = para.kvx;
wa  = para.wa;

A = [[ -(cos(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta), (sin(phi)*(b_gy - gyro_y))/(sin(theta)^2 - 1),    0,    0,  -1, -(sin(phi)*sin(theta))/cos(theta),   0,   0]
     [                          sin(phi)*(b_gy - gyro_y),                                             0,    0,    0,   0,                         -cos(phi),   0,   0]
     [                                                 0,                                  g*cos(theta), -kvx,    0,   0,                                 0,   0,   0]
     [                            -g*cos(phi)*cos(theta),                         g*sin(phi)*sin(theta),    0, -kvy,   0,                                 0,   0,   0]
     [                                                 0,                                             0,    0,    0, -wg,                                 0,   0,   0]
     [                                                 0,                                             0,    0,    0,   0,                               -wg,   0,   0]
     [                                                 0,                                             0,    0,    0,   0,                                 0, -wa,   0]
     [                                                 0,                                             0,    0,    0,   0,                                 0,   0, -wa]];

% [-(cos(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta), (sin(phi)*(b_gy - gyro_y))/(sin(theta)^2 - 1),    0,    0,  -1, -(sin(phi)*sin(theta))/cos(theta),   0,   0]
% [                         sin(phi)*(b_gy - gyro_y),                                             0,    0,    0,   0,                         -cos(phi),   0,   0]
% [                                                0,                                  g*cos(theta), -kvx,    0,   0,                                 0,   0,   0]
% [                           -g*cos(phi)*cos(theta),                         g*sin(phi)*sin(theta),    0, -kvy,   0,                                 0,   0,   0]
% [                                                0,                                             0,    0,    0, -wg,                                 0,   0,   0]
% [                                                0,                                             0,    0,    0,   0,                               -wg,   0,   0]
% [                                                0,                                             0,    0,    0,   0,                                 0, -wa,   0]
% [                                                0,                                             0,    0,    0,   0,                                 0,   0, -wa]

F_k = eye(8) + para.Ts*A;

return

function H_k = calc_H_k(x, u, para)

phi   = x(1);
theta = x(2);
vx    = x(3);
vy    = x(4);
b_gx  = x(5);
b_gy  = x(6);
b_ax  = x(7);
b_ay  = x(8);

gyro_x = u(1);
gyro_y = u(2);

wg  = para.wg;
g   = para.g;
kvx = para.kvx;
kvy = para.kvx;
wa  = para.wa;

H_k = [[ 0, 0, -kvx,    0, 0, 0, 1, 0]
       [ 0, 0,    0, -kvy, 0, 0, 0, 1]];

% [0, 0, -kvx,    0, 0, 0, 1, 0]
% [0, 0,    0, -kvy, 0, 0, 0, 1]
       
return

function [Qk, Rk] = calc_Q_R_ekf_test(x, u, para, var_fx, var_gy, rho)

phi   = x(1);
theta = x(2);
vx    = x(3);
vy    = x(4);
b_gx  = x(5);
b_gy  = x(6);
b_ax  = x(7);
b_ay  = x(8);

gyro_x = u(1);
gyro_y = u(2);

wg  = para.wg;
g   = para.g;
kvx = para.kvx;
kvy = para.kvx;
wa  = para.wa;

var_fx_0 = var_fx(1);
var_fx_1 = var_fx(2);
var_fx_2 = var_fx(3);
var_fx_3 = var_fx(4);
var_fx_4 = var_fx(5);
var_fx_5 = var_fx(6);
var_fx_6 = var_fx(7);
var_fx_7 = var_fx(8);

var_gy_0 = var_gy(1);
var_gy_1 = var_gy(2);

Q = [[ var_fx_0 + (var_fx_1*sin(phi)^2*sin(theta)^2)/cos(theta)^2, (var_fx_1*cos(phi)*sin(phi)*sin(theta))/cos(theta),        0,        0,        0,        0,        0,        0]
     [         (var_fx_1*cos(phi)*sin(phi)*sin(theta))/cos(theta),                                var_fx_1*cos(phi)^2,        0,        0,        0,        0,        0,        0]
     [                                                          0,                                                  0, var_fx_2,        0,        0,        0,        0,        0]
     [                                                          0,                                                  0,        0, var_fx_3,        0,        0,        0,        0]
     [                                                          0,                                                  0,        0,        0, var_fx_4,        0,        0,        0]
     [                                                          0,                                                  0,        0,        0,        0, var_fx_5,        0,        0]
     [                                                          0,                                                  0,        0,        0,        0,        0, var_fx_6,        0]
     [                                                          0,                                                  0,        0,        0,        0,        0,        0, var_fx_7]];

% [var_fx_0 + (var_fx_1*sin(phi)^2*sin(theta)^2)/cos(theta)^2, (var_fx_1*cos(phi)*sin(phi)*sin(theta))/cos(theta),        0,        0,        0,        0,        0,        0]
% [        (var_fx_1*cos(phi)*sin(phi)*sin(theta))/cos(theta),                                var_fx_1*cos(phi)^2,        0,        0,        0,        0,        0,        0]
% [                                                         0,                                                  0, var_fx_2,        0,        0,        0,        0,        0]
% [                                                         0,                                                  0,        0, var_fx_3,        0,        0,        0,        0]
% [                                                         0,                                                  0,        0,        0, var_fx_4,        0,        0,        0]
% [                                                         0,                                                  0,        0,        0,        0, var_fx_5,        0,        0]
% [                                                         0,                                                  0,        0,        0,        0,        0, var_fx_6,        0]
% [                                                         0,                                                  0,        0,        0,        0,        0,        0, var_fx_7]

% [~, A] = calc_F_k(x, u, para);

Qk = para.Ts * Q;
% Qk = para.Ts * ( Q + 1/2*( Q*A.' + A*Q )*para.Ts + 1/3*( A*Q*A.' )*para.Ts^2 );
% Qk = para.Ts * ( Q + 1/2*( Q*A.' + A*Q )*para.Ts + 1/3*( 1/2*(Q*A.'*A.' + A*A*Q) + A*Q*A.' )*para.Ts^2 + 1/8*(A*Q*A.'*A.' + A*A*Q*A.')*para.Ts^3 + 1/20*A*A*Q*A.'*A.'*para.Ts^4);
 
Rk = 1/para.Ts*[[ rho*var_gy_0,            0]
                [            0, rho*var_gy_1]];

% [rho*var_gy_0,            0]
% [           0, rho*var_gy_1]
         
return

function x_k = fxd(x, u, para)

phi   = x(1);
theta = x(2);
vx    = x(3);
vy    = x(4);
b_gx  = x(5);
b_gy  = x(6);
b_ax  = x(7);
b_ay  = x(8);

gyro_x = u(1);
gyro_y = u(2);

wg  = para.wg;
g   = para.g;
kvx = para.kvx;
kvy = para.kvx;
wa  = para.wa;

x_k = x + para.Ts*[ gyro_x - b_gx - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)
                                                           -cos(phi)*(b_gy - gyro_y)
                                                               g*sin(theta) - kvx*vx
                                                    - kvy*vy - g*cos(theta)*sin(phi)
                                                                            -b_gx*wg
                                                                            -b_gy*wg
                                                                            -b_ax*wa
                                                                            -b_ay*wa];

% gyro_x - b_gx - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)
%                                        -cos(phi)*(b_gy - gyro_y)
%                                            g*sin(theta) - kvx*vx
%                                 - kvy*vy - g*cos(theta)*sin(phi)
%                                                         -b_gx*wg
%                                                         -b_gy*wg
%                                                         -b_ax*wa
%                                                         -b_ay*wa

return

function y_k = gy(x, u, para)

phi   = x(1);
theta = x(2);
vx    = x(3);
vy    = x(4);
b_gx  = x(5);
b_gy  = x(6);
b_ax  = x(7);
b_ay  = x(8);

gyro_x = u(1);
gyro_y = u(2);

wg  = para.wg;
g   = para.g;
kvx = para.kvx;
kvy = para.kvx;
wa  = para.wa;

y_k = [ b_ax - kvx*vx
        b_ay - kvy*vy];

% b_ax - kvx*vx
% b_ay - kvy*vy
    
return