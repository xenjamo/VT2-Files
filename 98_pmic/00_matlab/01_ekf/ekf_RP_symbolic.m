% %% assumes we measure -(k1x, k1y)/m*v in acc sensor, extend with acc bias (only acc sensor)

clc, clear all

syms gyro_x gyro_y wg b_gx b_gy phi theta psi g kvx kvy vx vy wa

% kvx = kx/m
% kvy = ky/m
% wg = 1/tau_g

ang  = [phi; theta];
v    = [vx; vy];
b_g  = [b_gx; b_gy];
gyro = [gyro_x; gyro_y];
x    = [ang; v; b_g];

% gravitational field
grav = [0; 0; -g];

% earth to body rotation
CBE = [[                              cos(psi)*cos(theta),                              cos(theta)*sin(psi),         -sin(theta)]; ...
       [ cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(theta)*sin(phi)]; ...
       [ sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), cos(phi)*cos(theta)]];

% body to earth angle velocities
JEB = [[ 1, sin(phi)*sin(theta)/cos(theta), cos(phi)*sin(theta)/cos(theta)]; ...
       [ 0,                       cos(phi),                      -sin(phi)]; ...
       [ 0,            sin(phi)/cos(theta),            cos(phi)/cos(theta)]];

% syms gyro_z vz
% % test coriolis
% cross([gyro_x, gyro_y, gyro_z].' , [vx, vy, vz].')
% % gyro_y*vz - gyro_z*vy
% % gyro_z*vx - gyro_x*vz
% % gyro_x*vy - gyro_y*vx

% x_k = [ang; v; b_g; b_a]
fx = [JEB(1:2,1:2) * (gyro - b_g); ...
      CBE(1,:) * grav - kvx * vx; ...
      CBE(2,:) * grav - kvy * vy; ...
      -wg * b_g]

gy = [-kvx * vx; -kvy * vy]

Aa = simplify(jacobian(fx, x))
Ba = simplify(jacobian(fx, gyro))
Ca = simplify(jacobian(gy, x))
Da = simplify(jacobian(gy, gyro))

% [n_gyro; n_v; n_b_g; n_b_a]
Ga = [[-JEB(1:2,1:2), zeros(2,4)]; [zeros(4,2), eye(4,4)]];

syms var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5 var_fx_6 var_fx_7
syms var_gy_0 var_gy_1
var_fx = [var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5].';
var_gy = [var_gy_0 var_gy_1].';
syms Ts rho
Qa = simplify(Ga * diag(var_fx) * Ga.')
Ra = simplify( rho * diag(var_gy) )

% syms s1 c1 s2 c2 Ts x_0 x_1 x_2 x_3 x_4 x_5 u_0 u_1 x_6 x_7
% subs(Aa, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(eye(6) + Ts*Aa, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(Ca, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(Qa, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(x + Ts*fx, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.' v.' ang.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1 x_2 x_3 x_0 x_1])
% subs(gy, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.' v.' ang.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1 x_2 x_3 x_0 x_1])

Aa0 = subs(Aa, [x; gyro], zeros(8,1))
Ba0 = subs(Ba, [x; gyro], zeros(8,1))
Ca0 = subs(Ca, [x; gyro], zeros(8,1))
Da0 = subs(Da, [x; gyro], zeros(8,1))
Ga0 = subs(Ga, [x; gyro], zeros(8,1))

Ts = 1/50;
vars = [kvx kvy wg g];
num_vals = [0.2 0.2 2*pi*0.1*0 9.81]
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

% 6 states
% [n_gyro; n_v; n_b_g; n_b_a]
var_fx = [0.01*[1 1], 10*[1 1], 0.02*[1 1]];
% [n_acc]
var_gy = [5.0*[1 1]];
rho = 0.1;

Qn0 = Gn0 * diag(var_fx) * Gn0.'
;
% nx = size(An0,1);
% F = [[-An0, Qn0];[zeros(nx,nx), An0.']];
% Fd = expm(F*Ts);
% Qn0d = Fd(nx+1:end,nx+1:end).'*Fd(1:nx,nx+1:end)
Qn0d = Ts * Qn0;
% Qn0d = Ts * ( Qn0 + 1/2*( Qn0*An0.' + An0*Qn0 )*Ts + 1/3*( An0*Qn0*An0.' )*Ts^2 )
% Qn0d = Ts * ( Qn0 + 1/2*( Qn0*An0.' + An0*Qn0 )*Ts + 1/3*( 1/2*(Qn0*An0.'*An0.' + An0*An0*Qn0) + An0*Qn0*An0.' )*Ts^2 + ...
%               1/8*(An0*Qn0*An0.'*An0.' + An0*An0*Qn0*An0.')*Ts^3 + 1/20*An0*An0*Qn0*An0.'*An0.'*Ts^4)
Rn0d = rho * 1/Ts * diag(var_gy)

format long
K_k_stat = dlqr(An0d.', Cn0.', Qn0d, Rn0d).'
format short

eig_d = eig(An0d - K_k_stat*Cn0);
eig_c = log(eig_d)/Ts;
fn = abs(eig_c)/2/pi
Dn = cos(pi-angle(eig_c))

sys_d = ss(An0d - K_k_stat*Cn0, [-Bn0*Ts, -K_k_stat], eye(6), zeros(6,4), Ts);
% damp(sys_d)

figure(1)
sigma(sys_d), grid on, hold on
xlim([1e-5 1e2])
title('EKF Error Dynamics')

% figure(2)
% bodemag(sys_d([1 2 3 4],:)), grid on, hold on
% axis([1e-8 1e1 1e-6 1e2])
% title('EKF Error Dynamics')


%% assumes we measure -(k1x, k1y)/m*v in acc sensor, extend with acc bias (only acc sensor)

% % THIS IS THE ACTUAL C++ IMPLEMENTATION
% 
% clc, clear all
% 
% syms gyro_x gyro_y wg b_gx b_gy phi theta psi g kvx kvy vx vy wa b_ax b_ay
% 
% % kvx = kx/m
% % kvy = ky/m
% % wg = 1/tau_g
% % wa = 1/tau_a
% 
% ang  = [phi; theta];
% v    = [vx; vy];
% b_g  = [b_gx; b_gy];
% b_a  = [b_ax; b_ay];
% gyro = [gyro_x; gyro_y];
% x    = [ang; v; b_g; b_a];
% 
% % gravitational field
% grav = [0; 0; -g];
% 
% % earth to body rotation
% CBE = [[                              cos(psi)*cos(theta),                              cos(theta)*sin(psi),         -sin(theta)]; ...
%        [ cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(theta)*sin(phi)]; ...
%        [ sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), cos(phi)*cos(theta)]];
% 
% % body to earth angle velocities
% JEB = [[ 1, sin(phi)*sin(theta)/cos(theta), cos(phi)*sin(theta)/cos(theta)]; ...
%        [ 0,                       cos(phi),                      -sin(phi)]; ...
%        [ 0,            sin(phi)/cos(theta),            cos(phi)/cos(theta)]];
% 
% % x_k = [ang; v; b_g; b_a]
% fx = [JEB(1:2,1:2) * (gyro - b_g); ...
%       CBE(1,:) * grav - kvx * vx; ...
%       CBE(2,:) * grav - kvy * vy; ...
%       -wg * b_g; ...
%       -wa * b_a]
% 
% gy = [-kvx * vx + b_ax; -kvy * vy + b_ay]
% 
% Aa = simplify(jacobian(fx, x))
% Ba = simplify(jacobian(fx, gyro))
% Ca = simplify(jacobian(gy, x))
% Da = simplify(jacobian(gy, gyro))
% 
% % [n_gyro; n_v; n_b_g; n_b_a]
% Ga = [[-JEB(1:2,1:2), zeros(2,6)]; [zeros(6,2), eye(6,6)]];
% 
% syms var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5 var_fx_6 var_fx_7
% syms var_gy_0 var_gy_1
% var_fx = [var_fx_0 var_fx_1 var_fx_2 var_fx_3 var_fx_4 var_fx_5 var_fx_6 var_fx_7].';
% var_gy = [var_gy_0 var_gy_1].';
% syms Ts rho
% Qa = simplify(Ga * diag(var_fx) * Ga.')
% Ra = simplify( rho * diag(var_gy) )
% 
% syms s1 c1 s2 c2 Ts x_0 x_1 x_2 x_3 x_4 x_5 u_0 u_1 x_6 x_7
% subs(Aa, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(eye(8) + Ts*Aa, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(Ca, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(Qa, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1])
% subs(x + Ts*fx, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.' v.' ang.' b_a.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1 x_2 x_3 x_0 x_1  x_6 x_7])
% subs(gy, [sin(phi), cos(phi), sin(theta), cos(theta) b_g.' gyro.' v.' ang.' b_a.'], ...
%          [s1 c1 s2 c2 x_4 x_5 u_0 u_1 x_2 x_3 x_0 x_1 x_6 x_7])
% 
% Aa0 = subs(Aa, [x; gyro], zeros(10,1))
% Ba0 = subs(Ba, [x; gyro], zeros(10,1))
% Ca0 = subs(Ca, [x; gyro], zeros(10,1))
% Da0 = subs(Da, [x; gyro], zeros(10,1))
% Ga0 = subs(Ga, [x; gyro], zeros(10,1))
% 
% Ts = 1/50;
% vars = [kvx kvy wg g wa];
% num_vals = [0.25 0.25 2*pi*0.08 9.81 2*pi*0.16]
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
% % [n_gyro; n_v; n_b_g; n_b_a]
% var_fx = 0.008*[35.0*[1 1], 400*[1 1], 0.1*[1 1], 600*[1 1]];
% % [n_acc]
% var_gy = 8.0*[1 1];
% rho = 1.0;
% 
% Qn0 = Gn0 * diag(var_fx) * Gn0.';
% % nx = size(An0,1);
% % F = [[-An0, Qn0];[zeros(nx,nx), An0.']];
% % Fd = expm(F*Ts);
% % Qn0d = Fd(nx+1:end,nx+1:end).'*Fd(1:nx,nx+1:end)
% Qn0d = Ts * Qn0;
% % Qn0d = Ts * ( Qn0 + 1/2*( Qn0*An0.' + An0*Qn0 )*Ts + 1/3*( An0*Qn0*An0.' )*Ts^2 )
% % Qn0d = Ts * ( Qn0 + 1/2*( Qn0*An0.' + An0*Qn0 )*Ts + 1/3*( 1/2*(Qn0*An0.'*An0.' + An0*An0*Qn0) + An0*Qn0*An0.' )*Ts^2 + ...
% %               1/8*(An0*Qn0*An0.'*An0.' + An0*An0*Qn0*An0.')*Ts^3 + 1/20*An0*An0*Qn0*An0.'*An0.'*Ts^4)
% Rn0d = rho * 1/Ts * diag(var_gy)
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
% sys_d = ss(An0d - K_k_stat*Cn0, [-Bn0*Ts, -K_k_stat], eye(8), zeros(8,4), Ts);
% % damp(sys_d)
% 
% figure(1)
% sigma(sys_d), grid on, hold on
% xlim([1e-5 1e2])
% title('EKF Error Dynamics')
% 
% % figure(2)
% % bodemag(sys_d([1 2 3 4],:)), grid on, hold on
% % axis([1e-8 1e1 1e-6 1e2])
% % title('EKF Error Dynamics')
