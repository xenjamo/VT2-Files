clear all, clc; close all;

%%
load data.mat;

t = t - t(1);
lon = hpllh(:,1)*pi/180;
lat = hpllh(:,2)*pi/180;
h = hpllh(:,3);
figure(7)
plot3(lon , lat, h), grid on
title('llh'), xlabel('lon'), ylabel('lat'), zlabel('h')

%%

pos_ecef = transformWGS84ToECEF_R(lat, lon, h); % R stands for rad

figure(8)
plot3(pos_ecef(:,1) , pos_ecef(:,2) , pos_ecef(:,3)), grid on
axis equal, title('pos ecef'), xlabel('x'), ylabel('y'), zlabel('z')

%% HACK
ind0 = 4;
pos_ecef_0 = transformWGS84ToECEF_R(lat(ind0), lon(ind0), h(ind0));
phi = lon(ind0);
la = lat(ind0);
R_ecefToLocal_0 = [ -sin(phi),          cos(phi),       0; ...
                    -cos(phi)*sin(la), -sin(la)*sin(phi), cos(la); ...
                     cos(la)*cos(phi),  cos(la)*sin(phi), sin(la)];
pos_enu = (pos_ecef - pos_ecef_0) * R_ecefToLocal_0.';
                 
figure(9)
plot3(pos_enu(:,1) , pos_enu(:,2) , pos_enu(:,3)), grid on
axis equal, title('pos ecef t'), xlabel('x'), ylabel('y'), zlabel('z')