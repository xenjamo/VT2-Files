function [X_k, K_k, Q_k, R_k, P_k, F_k0, H_k0, P_svd_k, Q_k0, R_k0, P_k0, RPY_k, xi_k, Kmat_k] ...
           = ekf_RPYVP_17states_no_baro_vel_wrt_earth_delayed_01(Z, x_k, para, var_fxn, var_gxn, mag0, rho)

% acc for lin. model
Z0 = zeros(size(Z, 2));
Z0(6) = 9.81;

[Q_k, R_k] = calc_Q_R_ekf_test(zeros(13), zeros(3), para, var_fxn, var_gxn, rho);
F_k = calc_F_k(x_k, Z0, para);
H_k = calc_H_k(x_k, Z0, para, mag0);
P_k = 1 * para.scale_P0 * eye(17);

Q_k0 = Q_k;
R_k0 = R_k;
P_k0 = P_k;
F_k0 = F_k;
H_k0 = H_k;

X_k  = zeros(size(Z,1),length(x_k));
X_k(1,:) = x_k;
P_svd_k  = zeros(size(Z,1),17);
xi_k  = zeros(size(Z,1),1);
Kmat_k = zeros(size(Z,1),length(x_k)*7);

% delayed position estimates
nd_rspos = para.nd_rspos;
pos_past = zeros(nd_rspos+1, 3);
% pos_past = [px(k)   , py(k)   , pz(k)   ; ...
%             px(k-1) , py(k-1) , pz(k-1) ; ...
%                 .        .        .
%             px(k-nd), py(k-nd), pz(k-nd)]
   
for i = 1:size(X_k,1)
    
    F_k = calc_F_k(x_k, Z(i,:), para);
    H_k = calc_H_k(x_k, Z(i,:), para, mag0);
    [Q_k, R_k] = calc_Q_R_ekf_test(x_k, Z(i,:), para, var_fxn, var_gxn, rho);
    
    x_k = fxd(x_k, Z(i,:), para);
    P_k = F_k * P_k * F_k.' + Q_k;
        
    % update delayed position estimation
    pos_past = [x_k([7 8 9]).'; pos_past(1:nd_rspos,:)];
    
    x_k_past = x_k;
    x_k_past([7 8 9]) = pos_past(nd_rspos+1,:).';
    
    e = Z(i,[4 5 7 8 10 11 12]).' - gy(x_k_past, Z(i,:), para, mag0);
    
    S_inv_k = eye(7)/(H_k*P_k*H_k.' + R_k);
    K_k = (P_k*H_k.')*S_inv_k;
            
    x_k = x_k + K_k*e;
    P_k = (eye(17) - K_k*H_k) * P_k;
    
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
b_mx = x(16);
b_my = x(17);



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


acc_z  = u(6);

A = [[                                                                                                                                                                                                                                                                                                                                                                                                              (sin(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta) - (cos(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta),                                                                                                                                                                                                                                                                                                                                                                                                                                                        -(b_gz*cos(phi) - gyro_z*cos(phi) + b_gy*sin(phi) - gyro_y*sin(phi))/cos(theta)^2,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,     -1, -(sin(phi)*sin(theta))/cos(theta), -(cos(phi)*sin(theta))/cos(theta),      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                              cos(phi)*(b_gz - gyro_z) + sin(phi)*(b_gy - gyro_y),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                         -cos(phi),                          sin(phi),      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                    (sin(phi)*(b_gz - gyro_z))/cos(theta) - (cos(phi)*(b_gy - gyro_y))/cos(theta),                                                                                                                                                                                                                                                                                                                                                                                                                                - (cos(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta)^2 - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)^2,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,              -sin(phi)/cos(theta),              -cos(phi)/cos(theta),      0,      0,                                                  0,   0,   0]
[ (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(acc_z - b_az) - vy*(kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) + kvy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))) + 2*kvy*vx*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - kvy*vz*cos(theta)*(sin(psi) - 2*cos(phi)^2*sin(psi) + 2*cos(phi)*cos(psi)*sin(phi)*sin(theta)), vy*cos(theta)*(2*kvx*cos(psi)*sin(psi)*sin(theta) - kvy*cos(phi)*cos(psi)^2*sin(phi) + kvy*cos(phi)*sin(phi)*sin(psi)^2 - 2*kvy*cos(psi)*sin(phi)^2*sin(psi)*sin(theta)) - vz*(kvx*cos(psi)*sin(theta)^2 - kvx*cos(psi)*cos(theta)^2 + kvy*cos(psi)*cos(theta)^2*sin(phi)^2 - kvy*cos(psi)*sin(phi)^2*sin(theta)^2 + kvy*cos(phi)*sin(phi)*sin(psi)*sin(theta)) + cos(phi)*cos(psi)*cos(theta)*(acc_z - b_az) + 2*vx*cos(psi)*cos(theta)*(kvx*cos(psi)*sin(theta) + kvy*cos(phi)*sin(phi)*sin(psi) - kvy*cos(psi)*sin(phi)^2*sin(theta)), (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) - vz*(kvx*cos(theta)*sin(psi)*sin(theta) - kvy*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))) - vx*(2*kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - 2*kvx*cos(psi)*cos(theta)^2*sin(psi)) + vy*(kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))^2 - kvy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))^2 - kvx*cos(psi)^2*cos(theta)^2 + kvx*cos(theta)^2*sin(psi)^2),                                                       - kvy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))^2 - kvx*cos(psi)^2*cos(theta)^2, kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - kvx*cos(psi)*cos(theta)^2*sin(psi), kvx*cos(psi)*cos(theta)*sin(theta) + kvy*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)), 0, 0, 0,      0,                                 0,                                 0,      0,      0, - sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta),   0,   0]
[ 2*kvy*vy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - vx*(kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) + kvy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))) - (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(acc_z - b_az) - kvy*vz*cos(theta)*(2*cos(phi)^2*cos(psi) - cos(psi) + 2*cos(phi)*sin(phi)*sin(psi)*sin(theta)), vz*(kvx*cos(theta)^2*sin(psi) - kvx*sin(psi)*sin(theta)^2 - kvy*cos(theta)^2*sin(phi)^2*sin(psi) + kvy*sin(phi)^2*sin(psi)*sin(theta)^2 + kvy*cos(phi)*cos(psi)*sin(phi)*sin(theta)) + vx*cos(theta)*(2*kvx*cos(psi)*sin(psi)*sin(theta) - kvy*cos(phi)*cos(psi)^2*sin(phi) + kvy*cos(phi)*sin(phi)*sin(psi)^2 - 2*kvy*cos(psi)*sin(phi)^2*sin(psi)*sin(theta)) + cos(phi)*cos(theta)*sin(psi)*(acc_z - b_az) - 2*vy*cos(theta)*sin(psi)*(kvy*cos(phi)*cos(psi)*sin(phi) - kvx*sin(psi)*sin(theta) + kvy*sin(phi)^2*sin(psi)*sin(theta)), (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) + vz*(kvx*cos(psi)*cos(theta)*sin(theta) + kvy*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))) + vy*(2*kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - 2*kvx*cos(psi)*cos(theta)^2*sin(psi)) + vx*(kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))^2 - kvy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))^2 - kvx*cos(psi)^2*cos(theta)^2 + kvx*cos(theta)^2*sin(psi)^2), kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - kvx*cos(psi)*cos(theta)^2*sin(psi),                                                       - kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))^2 - kvx*cos(theta)^2*sin(psi)^2, kvx*cos(theta)*sin(psi)*sin(theta) - kvy*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)), 0, 0, 0,      0,                                 0,                                 0,      0,      0,   cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta),   0,   0]
[                                                                                                                                                                                                                         -cos(theta)*(acc_z*sin(phi) - b_az*sin(phi) - kvy*vy*cos(psi) + kvy*vx*sin(psi) + 2*kvy*vy*cos(phi)^2*cos(psi) - 2*kvy*vx*cos(phi)^2*sin(psi) + 2*kvy*vz*cos(phi)*cos(theta)*sin(phi) + 2*kvy*vx*cos(phi)*cos(psi)*sin(phi)*sin(theta) + 2*kvy*vy*cos(phi)*sin(phi)*sin(psi)*sin(theta)),                                                         vy*(kvx*cos(theta)^2*sin(psi) - kvx*sin(psi)*sin(theta)^2 + kvy*sin(phi)*sin(theta)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) - kvy*cos(theta)^2*sin(phi)^2*sin(psi)) - vx*(kvx*cos(psi)*sin(theta)^2 - kvx*cos(psi)*cos(theta)^2 + kvy*sin(phi)*sin(theta)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) + kvy*cos(psi)*cos(theta)^2*sin(phi)^2) - vz*(2*kvx*cos(theta)*sin(theta) - 2*kvy*cos(theta)*sin(phi)^2*sin(theta)) - cos(phi)*sin(theta)*(acc_z - b_az),                                                                                                                                                                                                                                                                                               vy*(kvx*cos(psi)*cos(theta)*sin(theta) + kvy*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))) - vx*(kvx*cos(theta)*sin(psi)*sin(theta) - kvy*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))),                                kvx*cos(psi)*cos(theta)*sin(theta) + kvy*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)),                                kvx*cos(theta)*sin(psi)*sin(theta) - kvy*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)),                                                                - kvy*cos(theta)^2*sin(phi)^2 - kvx*sin(theta)^2, 0, 0, 0,      0,                                 0,                                 0,      0,      0,                               -cos(phi)*cos(theta),   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              1,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                                 0,                                 0,      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              1,                                                                                                               0, 0, 0, 0,      0,                                 0,                                 0,      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               1, 0, 0, 0,      0,                                 0,                                 0,      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0, -wg_xy,                                 0,                                 0,      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                            -wg_xy,                                 0,      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                                 0,                             -wg_z,      0,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                                 0,                                 0, -wa_xy,      0,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                                 0,                                 0,      0, -wa_xy,                                                  0,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                                 0,                                 0,      0,      0,                                              -wa_z,   0,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                                 0,                                 0,      0,      0,                                                  0, -wm,   0]
[                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         0,                                                                                                                                              0,                                                                                                                                              0,                                                                                                               0, 0, 0, 0,      0,                                 0,                                 0,      0,      0,                                                  0,   0, -wm]];

F_k = eye(17) + para.Ts*A;

return

function H_k = calc_H_k(x, u, para, mag0)

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
b_mx = x(16);
b_my = x(17);




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


mx0 = mag0(1);
my0 = mag0(2);
mz0 = mag0(3);

H_k = [[                                                                                                                                                  0,           kvx*(vz*cos(theta) + vx*cos(psi)*sin(theta) + vy*sin(psi)*sin(theta)),                                                                           -kvx*cos(theta)*(vy*cos(psi) - vx*sin(psi)),                               -kvx*cos(psi)*cos(theta),                                -kvx*cos(theta)*sin(psi),           kvx*sin(theta), 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
[ kvy*vy*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) - kvy*vx*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) - kvy*vz*cos(phi)*cos(theta), -kvy*sin(phi)*(vx*cos(psi)*cos(theta) - vz*sin(theta) + vy*cos(theta)*sin(psi)), kvy*vx*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) + kvy*vy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)), kvy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)), -kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)), -kvy*cos(theta)*sin(phi), 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
[                                                                                                                                                  0,            - mz0*cos(theta) - mx0*cos(psi)*sin(theta) - my0*sin(psi)*sin(theta),                                                                              cos(theta)*(my0*cos(psi) - mx0*sin(psi)),                                                      0,                                                       0,                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
[          mx0*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)) - my0*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)) + mz0*cos(phi)*cos(theta),   sin(phi)*(mx0*cos(psi)*cos(theta) - mz0*sin(theta) + my0*cos(theta)*sin(psi)),     - mx0*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) - my0*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)),                                                      0,                                                       0,                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
[                                                                                                                                                  0,                                                                               0,                                                                                                                     0,                                                      0,                                                       0,                        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
[                                                                                                                                                  0,                                                                               0,                                                                                                                     0,                                                      0,                                                       0,                        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
[                                                                                                                                                  0,                                                                               0,                                                                                                                     0,                                                      0,                                                       0,                        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]];

return

function [Qk, Rk] = calc_Q_R_ekf_test(x, u, para, var_fx, var_gx, rho)

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
b_mx = x(16);
b_my = x(17);




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
var_fx_15 = var_fx(16);
var_fx_16 = var_fx(17);

var_gy_0 = var_gx(1);
var_gy_1 = var_gx(2);
var_gy_2 = var_gx(3);
var_gy_3 = var_gx(4);
var_gy_4 = var_gx(5);
var_gy_5 = var_gx(6);
var_gy_6 = var_gx(7);

Q = [[ (var_fx_0*cos(theta)^2 + var_fx_2*cos(phi)^2*sin(theta)^2 + var_fx_1*sin(phi)^2*sin(theta)^2)/cos(theta)^2, (cos(phi)*sin(phi)*sin(theta)*(var_fx_1 - var_fx_2))/cos(theta), -(sin(theta)*(var_fx_2 + var_fx_1*sin(phi)^2 - var_fx_2*sin(phi)^2))/(sin(theta)^2 - 1),        0,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                                            (cos(phi)*sin(phi)*sin(theta)*(var_fx_1 - var_fx_2))/cos(theta),            var_fx_1 - var_fx_1*sin(phi)^2 + var_fx_2*sin(phi)^2,                                       (sin(2*phi)*(var_fx_1 - var_fx_2))/(2*cos(theta)),        0,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                    -(sin(theta)*(var_fx_2 + var_fx_1*sin(phi)^2 - var_fx_2*sin(phi)^2))/(sin(theta)^2 - 1),               (sin(2*phi)*(var_fx_1 - var_fx_2))/(2*cos(theta)),              -(var_fx_2 + var_fx_1*sin(phi)^2 - var_fx_2*sin(phi)^2)/(sin(theta)^2 - 1),        0,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0, var_fx_3,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0, var_fx_4,        0,        0,        0,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0, var_fx_5,        0,        0,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0, var_fx_6,        0,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0, var_fx_7,        0,        0,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0, var_fx_8,        0,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0, var_fx_9,         0,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0, var_fx_10,         0,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0,         0, var_fx_11,         0,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0,         0,         0, var_fx_12,         0,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0,         0,         0,         0, var_fx_13,         0,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0, var_fx_14,         0,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0,         0, var_fx_15,         0]
[                                                                                                          0,                                                               0,                                                                                       0,        0,        0,        0,        0,        0,        0,        0,         0,         0,         0,         0,         0,         0, var_fx_16]];

% [~, A] = calc_F_k(x, u, para);

Qk = para.Ts * Q;
% Qk = para.Ts * ( Q + 1/2*( Q*A.' + A*Q )*para.Ts + 1/3*( A*Q*A.' )*para.Ts^2 );
% Qk = para.Ts * ( Q + 1/2*( Q*A.' + A*Q )*para.Ts + 1/3*( 1/2*(Q*A.'*A.' + A*A*Q) + A*Q*A.' )*para.Ts^2 + 1/8*(A*Q*A.'*A.' + A*A*Q*A.')*para.Ts^3 + 1/20*A*A*Q*A.'*A.'*para.Ts^4);
 
Rk = 1/para.Ts*[[ rho*var_gy_0,            0,            0,            0,            0,            0,            0]
[            0, rho*var_gy_1,            0,            0,            0,            0,            0]
[            0,            0, rho*var_gy_2,            0,            0,            0,            0]
[            0,            0,            0, rho*var_gy_3,            0,            0,            0]
[            0,            0,            0,            0, rho*var_gy_4,            0,            0]
[            0,            0,            0,            0,            0, rho*var_gy_5,            0]
[            0,            0,            0,            0,            0,            0, rho*var_gy_6]];

return

function x_k = fxd(x, u, para)

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
b_mx = x(16);
b_my = x(17);




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


acc_z  = u(6);

x_k = x + para.Ts*[                                                                                                                                                                                                                                                                                                                          gyro_x - b_gx - (cos(phi)*sin(theta)*(b_gz - gyro_z))/cos(theta) - (sin(phi)*sin(theta)*(b_gy - gyro_y))/cos(theta)
                                                                                                                                                                                                                                                                                                                                                                                          sin(phi)*(b_gz - gyro_z) - cos(phi)*(b_gy - gyro_y)
                                                                                                                                                                                                                                                                                                                                                              - (sin(phi)*(b_gy - gyro_y))/cos(theta) - (cos(phi)*(b_gz - gyro_z))/cos(theta)
 (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(acc_z - b_az) - vx*(kvy*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))^2 + kvx*cos(psi)^2*cos(theta)^2) + vz*(kvx*cos(psi)*cos(theta)*sin(theta) + kvy*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))) + vy*(kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - kvx*cos(psi)*cos(theta)^2*sin(psi))
 vz*(kvx*cos(theta)*sin(psi)*sin(theta) - kvy*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))) - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(acc_z - b_az) - vy*(kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))^2 + kvx*cos(theta)^2*sin(psi)^2) + vx*(kvy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - kvx*cos(psi)*cos(theta)^2*sin(psi))
                                                                                                   vx*(kvx*cos(psi)*cos(theta)*sin(theta) + kvy*cos(theta)*sin(phi)*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))) - vz*(kvy*cos(theta)^2*sin(phi)^2 + kvx*sin(theta)^2) - g + vy*(kvx*cos(theta)*sin(psi)*sin(theta) - kvy*cos(theta)*sin(phi)*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))) + cos(phi)*cos(theta)*(acc_z - b_az)
                                                                                                                                                                                                                                                                                                                                                                                                                                           vx
                                                                                                                                                                                                                                                                                                                                                                                                                                           vy
                                                                                                                                                                                                                                                                                                                                                                                                                                           vz
                                                                                                                                                                                                                                                                                                                                                                                                                                  -b_gx*wg_xy
                                                                                                                                                                                                                                                                                                                                                                                                                                  -b_gy*wg_xy
                                                                                                                                                                                                                                                                                                                                                                                                                                   -b_gz*wg_z
                                                                                                                                                                                                                                                                                                                                                                                                                                  -b_ax*wa_xy
                                                                                                                                                                                                                                                                                                                                                                                                                                  -b_ay*wa_xy
                                                                                                                                                                                                                                                                                                                                                                                                                                   -b_az*wa_z
                                                                                                                                                                                                                                                                                                                                                                                                                                     -b_mx*wm
                                                                                                                                                                                                                                                                                                                                                                                                                                     -b_my*wm];
return

function y_k = gy(x, u, para, mag0)

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
b_mx = x(16);
b_my = x(17);




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


acc_z  = u(6);

mx0 = mag0(1);
my0 = mag0(2);
mz0 = mag0(3);

y_k = [                                                                        b_ax + kvx*vz*sin(theta) - kvx*vx*cos(psi)*cos(theta) - kvx*vy*cos(theta)*sin(psi)
 b_ay + kvy*vx*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) - kvy*vy*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) - kvy*vz*cos(theta)*sin(phi)
                                                                                 b_mx - mz0*sin(theta) + mx0*cos(psi)*cos(theta) + my0*cos(theta)*sin(psi)
          b_my - mx0*(cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta)) + my0*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta)) + mz0*cos(theta)*sin(phi)
                                                                                                                                                        px
                                                                                                                                                        py
                                                                                                                                                        pz];
                                                                                                                     
return