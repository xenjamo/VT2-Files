%%

clc, clear variables
%% read file
fid = fopen('bin_files/005.bin');
M = []; % initalize array
A = [];
A2 = [];
B = [];
C = [];
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
    % M = [M;M_];                     % write array to matrix
    A = [A;A_];
    A2 = [A2;A2_];
    B = [B;B_];
    C = [C;C_];
    D = [D;bitget(D_,8:-1:1)];
end
M = [A,A2,B,C,D];

fclose(fid);                        % close file

%% extract relevant data
size(M)

n = B(:,2);
t = B(:,1)/1000000;
dt = [0;diff(t)];
dn = [0;diff(n)];
t = t-t(1);
t0 = 0;
tend = t(end);

itow = B(:,3);
ditow = [0;diff(itow)];
ttff = B(1,4);
msss = B(:,5);
dmsss = [0;diff(msss)];
rtk_flags = B(:,6);

gnss_fix = C(:,1);
fix_type = C(:,2);
numSV = C(:,3); % number of visible satellites
lastCorrectionAge = C(:,4);


gyro = A(:,1:3);
accel = A(:,4:6);
mag = A(:,7:9);

% A_mag =  [0.9845097,  0.0000000,  0.0000000;
%           0.0127330,  0.9886070,  0.0000000;
%           0.0528879,  0.0021035,  1.0268833];
% b_mag =  [0.4526199, -0.3468273,  0.4043421];
% mag_calib = (mag - b_mag)*A_mag.';

quat = A(:,10:13);
rpy = A(:,14:16);

posCovNN = A(:,18);
posCovNE = A(:,19);
posCovND = A(:,20);
posCovEE = A(:,21);
posCovED = A(:,22);
posCovDD = A(:,23);
velCovNN = A(:,24);
velCovNE = A(:,25);
velCovND = A(:,26);
velCovEE = A(:,27);
velCovED = A(:,28);
velCovDD = A(:,29);

K_pos = [posCovNN, posCovNE, posCovND, posCovEE, posCovED, posCovDD];

K_vel = [velCovNN, velCovNE, velCovND, velCovEE, velCovED, velCovDD];

DOP = A(:,30:36);

relposNED = [A(:,37), A(:,38), A(:,39)];
relaccNED = [A(:,40), A(:,41), A(:,42)];

hMSL = A(:,43);
hAcc = A(:,44);
vAcc = A(:,45);

% llh = [A(:,46), A(:,47), A(:,48)];
relvelNED = [A(:,46), A(:,47), A(:,48)];

sAcc = A(:,49);
gSpeed = A(:,50);
headMot = A(:,51)*pi/180;
headAcc = A(:,52)*pi/180;
magDec = A(:,53)*pi/180;
magAcc = A(:,54)*pi/180;

hpllh = A2(:,1:3);                       % close file
%%

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
