clc, clear variables
%%

syms phi theta psi vx vy vz px py pz b_az b_mx b_my
syms gyro_x gyro_y gyro_z mx0 my0 acc_x acc_y acc_z
syms kvx kvy g
syms wm wa_z

x = [phi; theta; psi; vx; vy; vz; px; py; pz; b_az; b_mx; b_my];
b_m  = [b_mx; b_my];
gyro = [gyro_x; gyro_y; gyro_z];

% gravitational field
grav = [0; 0; -g];

% initial magnetic field vector (at startup of copter)
m0   = [mx0; my0];

% earth to body rotation
CBE = [[                              cos(psi)*cos(theta),                              cos(theta)*sin(psi),         -sin(theta)]; ...
       [ cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(theta)*sin(phi)]; ...
       [ sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), cos(phi)*cos(theta)]];
   
CEB = CBE.';

% body to earth angle velocities
JEB = [[ 1, sin(phi)*sin(theta)/cos(theta), cos(phi)*sin(theta)/cos(theta)]; ...
       [ 0,                       cos(phi),                      -sin(phi)]; ...
       [ 0,             sin(phi)/cos(theta),           cos(phi)/cos(theta)]];

% use indempotent matrix, see: Comparison of Attitude Estimation Techniques for Low-cost Unmanned Aerial Vehicles
mu = [0 0 1].'
Cmu = eye(3) - mu*mu.'
CBEmu = simplify(CBE*Cmu)
CBEmuy = simplify(CBE*Cmu*CBE.') % simplify(CBE*Cmu*CEB)
   
% % dvx_c/dt = vy*wz - vz*wy
% % dvy_c/dt = vz*wx - vx*wz
% % dvz_c/dt = vx*wy - vy*wx
% coriolis = [ vy * (gyro_z - b_gz) - vz * (gyro_y - b_gy); ...
%              vz * (gyro_x - b_gx) - vx * (gyro_z - b_gz);
%              vx * (gyro_y - b_gy) - vy * (gyro_x - b_gx)] % bias also in coriolis term
% % coriolis = [ vy * (gyro_z) - vz * (gyro_y); ...
% %              vz * (gyro_x) - vx * (gyro_z);
% %              vx * (gyro_y) - vy * (gyro_x)] % bias not in coriolis term

% x = [phi; theta; Bvx; Bvy; b_wx; b_wy; ...    | input: gyro_x, gyro_y, gyro_z | output: acc_x, acc_y
%      psi; b_wz ; b_mx; b_my; ...              | input: no additional input    | output: mag_x, mag_y, mag_z
%      Epz; Evz; b_az; ...                      | input: CEB(3,:)   * acc       | output: pz_gnss, ( baro )
%      Epx; Epy];                               | input: CEB(1:2,:) * acc
fx = [JEB * gyro; ...
      CBE * grav - diag([kvx; kvy; 0]) * [vx; vy; vz] + [0; 0; acc_z - b_az];
      CEB * [vx; vy; vz]; ...
      -wa_z * b_az; ...
      -wm * b_m]
  
gy = [-kvx * vx; ...
      -kvy * vy; ...
      CBEmu(1:2,1:2)*(m0 + b_m); ...
      [px, py, pz].']

Aa = simplify(jacobian(fx, x))
Ba = simplify(jacobian(fx, [gyro; acc_z]))
Ca = simplify(jacobian(gy, x))
Da = simplify(jacobian(gy, [gyro; acc_z]))

% [n_gyro; n_v; n_acc_z; n_pxy; n_pz; n_b_az; n_b_m]
Ga = [[-JEB, zeros(3,12-3)]; [zeros(12-3,3), eye(12-3,12-3)]];

% [n_acc; n_mag; n_pxy; n_pz]
Ua(3:4,3:5) = CBEmuy(1:2,:);
Ua(1:2,1:2) = eye(2);
Ua(5:7,6:8) = eye(3);

syms var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5 var_fx_6 var_fx_7 var_fx_8 var_fx_9 var_fx_10 var_fx_11
syms var_gy_0 var_gy_1 var_gy_2 var_gy_3 var_gy_4 var_gy_5 var_gy_6 var_gy_7
var_fx = [var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5 var_fx_6 var_fx_7 var_fx_8 var_fx_9 var_fx_10 var_fx_11]
var_gy = [var_gy_0 var_gy_1 var_gy_2 var_gy_3 var_gy_4 var_gy_5 var_gy_6 var_gy_7].';
syms Ts rho
Qa = simplify(Ga * diag(var_fx) * Ga.')
Ra = rho * simplify(Ua * diag(var_gy) * Ua.')

psi0 = 0*pi/180;
Aa0 = subs(Aa, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
Ba0 = subs(Ba, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
Ca0 = subs(Ca, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
Da0 = subs(Da, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
Ga0 = subs(Ga, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])

Ts = 1/50;
mag0 = [0.023271637037396  -0.166823804378510   0.453394383192062];
vars = [kvx kvy g wm wa_z mx0 my0 acc_z];
num_vals = [0.2 0.2 9.81 2*pi*0.8 0 mag0(1:2) 9.81]
An0 = double(subs(Aa0, vars, num_vals))
Bn0 = double(subs(Ba0, vars, num_vals))
Cn0 = double(subs(Ca0, vars, num_vals))
Dn0 = double(subs(Da0, vars, num_vals))
Gn0 = double(subs(Ga0, vars, num_vals))

Ob = obsv(An0, Cn0);
rank(Ob)
cond(Ob)

An0d = Ts*An0 + eye(size(An0));
Obd = obsv(An0d, Cn0);
rank(Obd)
cond(Obd)

eig(An0)

% [n_gyro; n_v; n_acc_z; n_pxy; n_pz; n_b_az; n_b_m]
var_fxn = [0.01*[1 1 1], 10*[1 1], 20, 1*[1 1], 1, 10, 4*[1 1]];

% [n_accxy; n_mag; n_pos]
var_gyn = [5.0*[1 1], 1.0*[1 1], 3*[1 1 1]];
rho = 0.1;

Qn0 = Gn0 * diag(var_fxn) * Gn0.';
Qn0d = Ts * Qn0;
Rn0d = rho * 1/Ts * diag(var_gyn)

format long
K_k_stat = dlqr(An0d.', Cn0.', Qn0d, Rn0d).'
format short

eig_d = eig(An0d - K_k_stat*Cn0);
eig_c = log(eig_d)/Ts;
fn = abs(eig_c)/2/pi
Dn = cos(pi-angle(eig_c))

sys_d = ss(An0d - K_k_stat*Cn0, [-Bn0*Ts, -K_k_stat], eye(12), zeros(12, 11), Ts);
% damp(sys_d)

figure(1)
sigma(sys_d), grid on, hold on
xlim([1e-5 1e2])
title('EKF Error Dynamics')

% figure(2)
% bodemag(sys_d(1:5,:)), grid on, hold on
% axis([1e-8 1e1 1e-6 1e2])
% title('EKF Error Dynamics')

%% just a try, not working yet

% % pose and velocity
% syms phi theta psi         % roll, pitch, yaw 3-2-1 euler angles
% syms vx vy vz              % velocity         w.r.t. earth frame
% syms px py pz              % position         w.r.t. earth frame
% 
% % angular rate
% syms gyro_x gyro_y gyro_z  % gyro             w.r.t. body frame
% syms b_gx   b_gy   b_gz    % gyro bias        w.r.t. body frame
% 
% % acc measurement
% syms acc_x acc_y acc_z     % acc              w.r.t. body frame
% syms b_ax  b_ay  b_ay      % acc bias         w.r.t. body frame
% 
% % velocity measurement and wind speed
% syms vel_n vel_w vel_u     % velocity         w.r.t. earth frame (NWU)
% syms v_wx  v_wy  v_wz      % wind speed       w.r.t. earth frame (NWU)
% 
% % position measurement
% syms pos_n  pos_w  pos_u   % velocity         w.r.t. earth frame (NWU)
% 
% % magnetometer measurement
% syms mag_x  mag_y  mag_z   % mag              w.r.t. body frame
% syms mag0_x mag0_y         % mag when initialising, level and towards north
%                            % therefor body, earth and NEU frame are aligned
%                            % this is a wrong assumption!
% syms b_mx b_my             % mag bias         w.r.t. body frame
% 
% % states
% ang   = [phi theta psi].'
% v     = [vx vy vz].'
% p     = [px py pz].'
% v_w   = [v_wx v_wy v_wz].'
% b_g   = [b_gx b_gy b_gz].'
% b_a   = [b_ax b_ay b_ay].'
% b_m   = [b_mx b_my].'
% 
% x = [ang; v; p; v_w; b_g; b_a; b_m]
% 
% % meassurement at initialisation
% mag0    = [mag0_x mag0_y].'
% 
% % measurements
% gyro  = [gyro_x gyro_y gyro_z].'
% acc   = [acc_x acc_y acc_z].'
% mag   = [mag_x mag_y mag_z].'
% vel   = [vel_n  vel_w  vel_u].'
% pos   = [pos_n  pos_w  pos_u].'
% 
% % bias dynamics (tuning parameters)
% syms wg_xy wg_z wa_xy wa_z wm
% wg = [wg_xy wg_xy wg_z]
% wa = [wa_xy wa_xy wa_z]
% wm = [wm wm]
% 
% % gravitational field
% syms g
% grav = [0; 0; -g];
% 
% % initial magnetic field vector (at startup of copter)
% mag0 = [mag0_x; mag0_y];
% 
% % earth to body rotation
% CBE = [[                              cos(psi)*cos(theta),                              cos(theta)*sin(psi),         -sin(theta)]; ...
%        [ cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(theta)*sin(phi)]; ...
%        [ sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), cos(phi)*cos(theta)]];
%    
% CEB = CBE.';
% 
% % body to earth jacobian (angle velocities)
% JEB = [[ 1, sin(phi)*sin(theta)/cos(theta), cos(phi)*sin(theta)/cos(theta)]; ...
%        [ 0,                       cos(phi),                      -sin(phi)]; ...
%        [ 0,             sin(phi)/cos(theta),           cos(phi)/cos(theta)]];
% 
% % use indempotent matrix, see: Comparison of Attitude Estimation Techniques for Low-cost Unmanned Aerial Vehicles
% mu = [0 0 1].'
% Cmu = eye(3) - mu*mu.'
% CBEmu = simplify(CBE*Cmu)
% CBEmuy = simplify(CBE*Cmu*CBE.') % simplify(CBE*Cmu*CEB)
%    
% % x = [theta; v; p; v_w; b_g; b_a; b_m]
% fx = [JEB * (gyro - b_g); ...
%       CEB * (acc  - b_a);
%       v + v_w; ...
%       -diag(wg) * b_g; ...
%       -diag(wa) * b_a; ...
%       -diag(wm) * b_m]
%   
% gy = [-kvx * vx; ...
%       -kvy * vy; ...
%       CBEmu(1:2,1:2)*(mag0 + b_m); ...
%       [px, py, pz].']
% 
% Aa = simplify(jacobian(fx, x))
% Ba = simplify(jacobian(fx, [gyro; acc_z]))
% Ca = simplify(jacobian(gy, x))
% Da = simplify(jacobian(gy, [gyro; acc_z]))
% 
% % [n_gyro; n_v; n_acc_z; n_pxy; n_pz; n_b_az; n_b_m]
% Ga = [[-JEB, zeros(3,12-3)]; [zeros(12-3,3), eye(12-3,12-3)]];
% 
% % [n_acc; n_mag; n_pxy; n_pz]
% Ua(3:4,3:5) = CBEmuy(1:2,:);
% Ua(1:2,1:2) = eye(2);
% Ua(5:7,6:8) = eye(3);
% 
% syms var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5 var_fx_6 var_fx_7 var_fx_8 var_fx_9 var_fx_10 var_fx_11
% syms var_gy_0 var_gy_1 var_gy_2 var_gy_3 var_gy_4 var_gy_5 var_gy_6 var_gy_7
% var_fx = [var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5 var_fx_6 var_fx_7 var_fx_8 var_fx_9 var_fx_10 var_fx_11]
% var_gy = [var_gy_0 var_gy_1 var_gy_2 var_gy_3 var_gy_4 var_gy_5 var_gy_6 var_gy_7].';
% syms Ts rho
% Qa = simplify(Ga * diag(var_fx) * Ga.')
% Ra = rho * simplify(Ua * diag(var_gy) * Ua.')
% 
% psi0 = 0*pi/180;
% Aa0 = subs(Aa, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
% Ba0 = subs(Ba, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
% Ca0 = subs(Ca, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
% Da0 = subs(Da, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
% Ga0 = subs(Ga, [x; gyro], [zeros(2,1); psi0; zeros(12,1)])
% 
% Ts = 1/50;
% mag0 = [0.023271637037396  -0.166823804378510   0.453394383192062];
% vars = [kvx kvy g wm wa_z mx0 my0 acc_z];
% num_vals = [0.2 0.2 9.81 2*pi*0.8 0 mag0(1:2) 9.81]
% An0 = double(subs(Aa0, vars, num_vals))
% Bn0 = double(subs(Ba0, vars, num_vals))
% Cn0 = double(subs(Ca0, vars, num_vals))
% Dn0 = double(subs(Da0, vars, num_vals))
% Gn0 = double(subs(Ga0, vars, num_vals))
% 
% Ob = obsv(An0, Cn0);
% rank(Ob)
% cond(Ob)
% 
% An0d = Ts*An0 + eye(size(An0));
% Obd = obsv(An0d, Cn0);
% rank(Obd)
% cond(Obd)
% 
% eig(An0)
% 
% % [n_gyro; n_v; n_acc_z; n_pxy; n_pz; n_b_az; n_b_m]
% var_fxn = [0.01*[1 1 1], 10*[1 1], 20, 1*[1 1], 1, 10, 4*[1 1]];
% 
% % [n_accxy; n_mag; n_pos]
% var_gyn = [5.0*[1 1], 1.0*[1 1], 3*[1 1 1]];
% rho = 0.1;
% 
% Qn0 = Gn0 * diag(var_fxn) * Gn0.';
% Qn0d = Ts * Qn0;
% Rn0d = rho * 1/Ts * diag(var_gyn)
% 
% format long
% K_k_stat = dlqr(An0d.', Cn0.', Qn0d, Rn0d).'
% format short
% 
% eig_d = eig(An0d - K_k_stat*Cn0);
% eig_c = log(eig_d)/Ts;
% fn = abs(eig_c)/2/pi
% Dn = cos(pi-angle(eig_c))
% 
% sys_d = ss(An0d - K_k_stat*Cn0, [-Bn0*Ts, -K_k_stat], eye(12), zeros(12, 11), Ts);
% % damp(sys_d)
% 
% figure(1)
% sigma(sys_d), grid on, hold on
% xlim([1e-5 1e2])
% title('EKF Error Dynamics')
% 
% % figure(2)
% % bodemag(sys_d(1:5,:)), grid on, hold on
% % axis([1e-8 1e1 1e-6 1e2])
% % title('EKF Error Dynamics')