clear all, close all, clc
%%
load data.mat

%%
dt = [0;diff(t)];
dn = [0;diff(n)];
n = n-n(1);
t = t-t(1);
t0 = 0;
tend = t(end);

figure(8);
subplot(311);
plot(n, relposNED(:,1))
title("relposN");
%xlim([t0 tend])
grid on;
subplot(312);
plot(n, relposNED(:,2))
title("relposE");
%xlim([t0 tend])
grid on;
subplot(313);
plot(n, relposNED(:,3))
title("relposD");
%xlim([t0 tend])
grid on;

figure(2);
plot(n,accel);
title("accel")
grid on;

%% 

p1_ind = 20:670;
p2_ind = 1267:1957;
p3_ind = 2836:3200;
p4_ind = 4450:5020;

p1 = relposNED(p1_ind,:);
p2 = relposNED(p2_ind,:);
p3 = relposNED(p3_ind,:);
p4 = relposNED(p4_ind,:);

p1_avg = [mean(p1(:,1)),mean(p1(:,2))];
p2_avg = [mean(p2(:,1)),mean(p2(:,2))];
p3_avg = [mean(p3(:,1)),mean(p3(:,2))];
p4_avg = [mean(p4(:,1)),mean(p4(:,2))];

l1 = p1_avg-p2_avg;
l2 = p4_avg-p3_avg;

v1 = [l1(2),l1(1)];
v1 = v1/norm(v1);
phi = atan2(v1(2),v1(1))*180/pi
v2 = [l2(2),l2(1)];
v2 = v2/norm(v2);
gamma = atan2(v2(2),v2(1))*180/pi

%%
figure(80)
subplot(221)
plot(p1(:,2), p1(:,1), '.');
hold on;
plot(p1_avg(2), p1_avg(1),'+');
grid on; hold off;
title("p1");
subplot(222)
plot(p2(:,2), p2(:,1), '.');
hold on;
plot(p2_avg(2), p2_avg(1),'+');
grid on; hold off;
title("p2");
subplot(223)
plot(p3(:,2), p3(:,1), '.');
hold on;
plot(p3_avg(2), p3_avg(1),'+');
grid on; hold off;
title("p3");
subplot(224)
plot(p4(:,2), p4(:,1), '.');
hold on;
plot(p4_avg(2), p4_avg(1),'+');
grid on; hold off;
title("p4");

figure(81);
plot([0,v1(1)], [0,v1(2)])
hold on; grid on;
plot([0,v2(1)], [0,v2(2)])
hold off;
xlim([-1 1])
ylim([-1 1])
title("vector")


