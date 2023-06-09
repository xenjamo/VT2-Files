clear all; clc, close all;
%% from ref.m
% these values dont add to 90 deg
% angle from x axis (east axis)
phi = -22.060361206351129;
gamma = 68.015863207738889;
% overwrite phi to get orthogonal grid
phi = gamma - 90;

%%
load data.mat
n = n-n(1);
% headMot is angle rel to North so we offset that by 90deg
% headMot = headMot + 90*pi/180; % no we dont
myWrapToPi = @(ang) atan2( sin(ang), cos(ang) );
headMot = myWrapToPi(headMot);
%%
figure(1)
plot(n,[headMot,headAcc])

grid on;

h1_ind = 504:834;
h2_ind = 970:1024;

h1 = mean(headMot(h1_ind)) * 180/pi;
h2 = mean(headMot(h2_ind)) * 180/pi;

%% mahony heading

A_mag =  [0.9845097,  0.0000000,  0.0000000;
          0.0127330,  0.9886070,  0.0000000;
          0.0528879,  0.0021035,  1.0268833];
b_mag =  [0.4526199, -0.3468273,  0.4043421];
mag_calib = (mag - b_mag)*A_mag.';

data = [gyro, accel, mag_calib, t, quat, rpy, headMot];

time = data(:,10);
time = time - time(1);
Ts = median(diff(time));
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
ind_quat = 11:14;
ind_rpy  = 15:17;
ind_headMot = 18;

% bessel
p = 2;         % pole at p rad/s
kp = 2 * p;
ki = kp^2 / 3;

para.kp = kp;
para.ki = ki;

rpy0 = 0 * [60, -60, 0] * 180/pi;
quat0 = rpy2quat(rpy0).';

[quatRP , biasRP ] = mahonyRP (data(:,ind_gyro), data(:,ind_acc), para, Ts, quat0);
[quatRPY, biasRPY] = mahonyRPY(data(:,ind_gyro), data(:,ind_acc), data(:,ind_mag), para, Ts, quat0);

rpyRP  = quat2rpy(quatRP );
rpyRPY = quat2rpy(quatRPY);
%%

figure(2)
plot(time, -headMot *180/pi)
hold on; grid on;
plot(time, rpyRPY(:,3) *180/pi)
hold off;
title("heading");