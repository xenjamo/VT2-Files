clc, clear all
addpath ..\99_fcn_bib\
%%

ex = [1, 0, 0].'
ez = [0, 0, 1].';

N = 1e2 + 1;
ang   = (0:N-1).' * 2*pi / (N-1) - pi;
theta_dot   = zeros(N, 1);
theta_cross = zeros(N, 1);
for i = 1:N
    v = calcRn(ez, ang(i)) * ex;
    theta_dot(i)   = acos(         dot(v, ex) );
    theta_cross(i) = asin( norm( cross(v, ex) ) );
end

figure(1)
plot(ang, [theta_dot, theta_cross]), grid on

