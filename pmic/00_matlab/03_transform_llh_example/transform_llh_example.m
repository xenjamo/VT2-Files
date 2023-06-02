load data_putty_05.mat

pos_ecef = transformWGS84ToECEF_R(data.lonlat(:,2), data.lonlat(:,1), data.height);

figure(8)
plot3(pos_ecef(:,1) , pos_ecef(:,2) , pos_ecef(:,3)), grid on, hold on
plot3(data.pos_ecef(:,1) , data.pos_ecef(:,2) , data.pos_ecef(:,3), 'color', [0 0.5 0]), grid on, hold off
axis equal, title('pos ecef'), xlabel('x'), ylabel('y'), zlabel('z')

% HACK
ind0 = 189;
pos_ecef_0 = transformWGS84ToECEF_R(data.lonlat(ind0,2), data.lonlat(ind0,1), data.height(ind0));
phi = data.lonlat(ind0,1);
la = data.lonlat(ind0,2);
R_ecefToLocal_0 = [ -sin(phi),          cos(phi),       0; ...
                    -cos(phi)*sin(la), -sin(la)*sin(phi), cos(la); ...
                     cos(la)*cos(phi),  cos(la)*sin(phi), sin(la)];
pos_enu = (pos_ecef - pos_ecef_0) * R_ecefToLocal_0.';
                 
figure(9)
plot(data.time(:), pos_enu(:,:)), grid on, hold on
plot(data.time(:), data.pos_enu(:,:)), grid on, ylim([-100 40]), hold off
xlabel('Time'), xlim([0 max(time)])