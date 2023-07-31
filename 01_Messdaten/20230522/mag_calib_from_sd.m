%%

clc, clear variables
%%
fid = fopen('mag_calib.bin');
M = []; % initalize array
A = [];
A2 = [];
B = [];
C_ = [];
D = [];

while 1
    A_ = fread(fid,[1 54],'float');
    A2_ = fread(fid,[1 3], 'double');
    B_ = fread(fid,[1 6],'uint32');
    C_ = fread(fid,[1 4],'uint8');
    D_ = fread(fid,[1 1],'uint8');
    M_ = [A_,A2_,B_,C_,D_];                   
    if size(M_) < 68                
        break;                      
    end
    M = [M;M_];                     % write array to matrix
    A = [A;A_];
    A2 = [A2;A2_];
    B = [B;B_];
    C_ = [C_;C_];
    D = [D;bitget(D_,8:-1:1)];
end

fclose(fid);                        % close file
%%
% data = readmatrix('putty_00.log'); % before calibration
% data = readmatrix('putty_01.log'); % after calibration
data = A(:,1:9);

% time = data(:,end) * 1e-3;
% time = time - time(1);
% data = data(:,1:end-1);
% Ts = median(diff(time));
% ind_gyro = 1:3;
% ind_acc  = 4:6;
ind_mag  = 7:9;

% use PROPER part of the measurement
mag = data(1:end,ind_mag);

C_ = [[0 0 1];[0 1 0];[1 0 0]];
[A_mag, b_mag] = MgnCalibration(mag*C_);
A_mag = C_.'*A_mag*C_;
expgfs0 = 1/mean(eig(A_mag)); % expected field strength
format long
A_mag = A_mag*expgfs0
b_mag = C_.'*b_mag

mag_calib = (mag - b_mag.')*A_mag.'; % m (Nx3) (N measpoints)

figure(1)
ax(1) = subplot(221);
plot(mag), grid on
title('uncalibrated')
ax(2) = subplot(223);
plot(sqrt(sum(mag.^2, 2))), grid on
ax(3) = subplot(222);
plot(mag_calib), grid on
title('calibrated')
ax(4) = subplot(224);
plot(sqrt(sum(mag_calib.^2, 2))), grid on
linkaxes(ax, 'x'), clear ax
xlim([0 size(mag, 1)])

figure(2)
plot3(mag(:,1),mag(:,2),mag(:,3), 'b.'), grid on, hold on
plot3(mag_calib(:,1),mag_calib(:,2),mag_calib(:,3), '.', 'color', [0 0.5 0]), hold off
axis equal
legend('uncalibrated', 'calibrated', 'location', 'best')

fprintf('A_mag << %10.7ff, %10.7ff, %10.7ff,\n', A_mag(1,:))
fprintf('         %10.7ff, %10.7ff, %10.7ff,\n', A_mag(2,:))
fprintf('         %10.7ff, %10.7ff, %10.7ff;\n', A_mag(3,:))
fprintf('b_mag << %10.7ff, %10.7ff, %10.7ff;\n', b_mag)
