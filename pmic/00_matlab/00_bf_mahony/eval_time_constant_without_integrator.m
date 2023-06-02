clc, clear variables

syms s kp ki
A = [0]
B = [1]
C = [1]
K = [kp]

eig(A - K*C)

Ggyro = simplify( C * (s - (A - K*C))^-1 * B )
Gacc  = simplify( C * (s - (A - K*C))^-1 * K )

%%
clc, clear all

p = 1;         % pole at p rad/s
kp = p;

s = tf('s');
Ggyro = 1/(kp + s);
Gacc  = kp/(kp + s);
pole(Ggyro)
damp(Ggyro)

figure(99)
bode(Ggyro * s, Gacc), hold on
bode(Ggyro * s + Gacc), hold off
