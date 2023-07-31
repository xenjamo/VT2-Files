clc, clear variables
%%

g = 9.81;
b_acc = 1;

var_gps = 0.001;
var_acc = 0.01;

Ts = 1e-3;
Tsim = 10 - Ts;
time = (0:Ts:Tsim).';

N = size(time, 1);

pos = step(tf(1, [1 1 1]), time);
vel = step(tf([1 0], [1 1 1]), time);
acc = step(tf([1 0 0], [1 1 1]), time);

y_gps = pos + sqrt(var_gps) * randn(N, 1);
y_acc = acc + b_acc + g + sqrt(var_acc) * randn(N, 1);

A = [0 1 0; ...
     0 0 -1; ...
     0 0 0];
B = [0; 1; 0];
C = [1 0 0];
A = eye(size(A)) + Ts * A;
B = Ts * B;

sys = ss(A, B, C, 0, Ts);

R = var_gps / Ts;
Q = B * B.' * var_acc * Ts;
Q(3,3) = 0.1 * Ts;
K = dlqr(A.', C.', Q, R).'; % statische Lösung Kalman-Filter
sys_est = ss(A - K*C, [B, K], eye(3), 0, Ts);

% eig_d = eig(sys_est.a);
% eig_c = log(eig_d) / Ts;
% fn = abs(eig_c) / 2 / pi
% Dn = cos(pi-angle(eig_c))
damp( sys_est )

x0_hat = [0, 0, b_acc].';
% x_hat_kf_static = lsim(sys_est, [y_acc - g, y_gps], time);
x_hat_kf_static = lsim(sys_est, [y_acc - g, y_gps], time, x0_hat);

%%

K_k = zeros(size(K));
P = 1 * eye(size(A));
x_k = x0_hat; % zeros(3,1);
x_hat_kf = zeros( size(x_hat_kf_static) );

for i = 1:N

    x_k = A * x_k + B * ( y_acc(i) - g );
    P = A * P * A.' + Q;
    
    e = y_gps(i) - C * x_k;
    K_k = (P * C.') / ( C * P * C.' + R );
    x_k = x_k + K_k * e;
    P = ( eye(size(A)) - K_k * C) * P;

    x_hat_kf(i,:) = x_k.';
end

figure(1)
subplot(311)
plot(time, [y_gps, x_hat_kf_static(:,1), x_hat_kf(:,1)]), grid on
subplot(312)
plot(time, [vel, x_hat_kf_static(:,2), x_hat_kf(:,2)]), grid on
subplot(313)
plot(time, y_acc - g), grid on % carefull here, just for visuallization

figure(2)
plot(time, [zeros(N,1), x_hat_kf_static(:,3), x_hat_kf(:,3)]), grid on

% figure(3)
% bode(sys_est(1,1), sys_est(1,2)), grid on


