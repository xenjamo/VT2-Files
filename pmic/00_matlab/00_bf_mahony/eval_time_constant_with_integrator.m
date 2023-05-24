clc, clear variables

syms s kp ki
A = [0 -1; 0 0]
B = [1; 0]
C = [1, 0]
K = [kp; -ki]

eig(A - K*C)

Ggyro = simplify( C * (s*eye(2) - (A - K*C))^-1 * B )
Gacc  = simplify( C * (s*eye(2) - (A - K*C))^-1 * K )

sol_real_poles = ...
    solve( kp^2 - 4*ki == 0, ki)
% kp^2/4

sol_compl_conj_poles = ...
    solve( (-kp^2 + 4*ki)^(1/2)/2 - kp/2 == 0, ki)
% kp^2/2

%%
clc, clear all

% Kp * (1 + 1/(Tn*s)) = kp + ki <->
% kp = Kp
% ki = Kp / Tn = kp / (1 / (2*pi*wn)) = kp * 2*pi*wn

% % just something
% kp = 2.0;
% ki = kp * (2*pi* 0.4);

% % double real pole
% p = 1;         % pole at p rad/s
% kp = 2 * p;
% ki = kp^2 / 4;

% % compl. conj. poles
% p = 1;         % pole at p rad/s
% kp = 2 * p;
% ki = kp^2 / 2;

% bessel
p = 1;         % pole at p rad/s
kp = 2 * p;
ki = kp^2 / 3;

s = tf('s');
Ggyro = s/(s^2 + kp*s + ki);
Gacc  = (ki + kp*s) / (s^2 + kp*s + ki);
pole(Ggyro)
damp(Ggyro)

figure(99)
bode(Ggyro * s, Gacc), hold on
bode(Ggyro * s + Gacc), hold off
