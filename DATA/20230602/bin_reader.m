clear all, clc, close all

%% read file
fid = fopen('001.bin');
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
t0 = 6;
tend = 303.5;

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
A_mag = [0.9858873,  0.0000000,  0.0000000;
         0.0105822,  0.9857083,  0.0000000;
         0.0572703,  0.0053311,  1.0284045];
b_mag = [0.4574184, -0.3294953,  0.3824048]';
mag_calib = (mag - b_mag.')*A_mag.';

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

relposNEU = [A(:,37), A(:,38), -A(:,39)];
relaccNEU = [A(:,40), A(:,41), A(:,42)];

hMSL = A(:,43);
hAcc = A(:,44);
vAcc = A(:,45);

llh = [A(:,46), A(:,47), A(:,48)];
relvelNEU = [A(:,49), A(:,50), -A(:,51)];

sAcc = A(:,52);
gSpeed = A(:,53);
headMot = A(:,54)*pi/180;

hpllh = A2(:,1:3);

save data.mat n t itow ttff msss rtk_flags...
     gnss_fix fix_type numSV lastCorrectionAge...
     gyro accel mag mag_calib quat rpy K_pos K_vel...
     DOP relposNEU relaccNEU hMSL hAcc vAcc llh...
     relvelNEU sAcc gSpeed headMot hpllh

%% plots

figure(999);
subplot(411)
plot(t,gnss_fix);
title("GPS Fix")
grid on;
subplot(412)
plot(t,fix_type);
title("Fix Type")
grid on;
subplot(413)
plot(t,numSV);
title("number of visible satellites")
ylim([0 max(numSV)+5])
grid on;
subplot(414)
plot(t,lastCorrectionAge);
title("Last RTK Correction Age")
grid on;

figure(9999);
subplot(411)
plot(t,dt);
title("dt Microcontroller")
grid on;
subplot(412)
plot(t,dn);
title("dSample")
grid on;
subplot(413)
plot(t,ditow);
title("dt GPS Time")
grid on;
subplot(414)
plot(t,dmsss);
title("dt RTC time")
grid on;

%%
figure(1);
plot(t,gyro);
title("gyro")
xlim([t0, tend])
grid on;

figure(2);
plot(t,accel);
title("accel")
xlim([t0, tend])
grid on;

figure(3);
subplot(211)
plot(t,mag);
title("mag")
xlim([t0, tend])
grid on;
subplot(212);
plot(t,mag_calib);
title("mag calib")
xlim([t0, tend])
grid on;

% carefull here
% Tmean = 6.0;
% ind0 = t > 1 & t < Tmean;
% mag0 = median(mag_calib(ind0)).';
% Ymag0 = atan2(mag_calib(ind0, 2), mag_calib(ind0, 1));
% Ymag0 = median(Ymag0);
% Ymag = atan2(mag_calib(:, 2), mag_calib(:,1));
% Ymag = atan2(sin(Ymag - Ymag0), cos(Ymag - Ymag0)); % unwrap(Ymag - Ymag0);
% figure(33)
% plot(t, unwrap([Ymag, heading])), grid on

figure(5);
plot(t,rpy);
title("rpy")
xlim([t0, tend])
grid on;
%%
figure(7);
plot3(relposNEU(:,2),relposNEU(:,1),relposNEU(:,3), '.-');
title("RELPOS")
grid on;

figure(8);
subplot(311);
plot(t, relposNEU(:,1))
title("relposN");
xlim([t0 tend])
grid on;
subplot(312);
plot(t, relposNEU(:,2))
title("relposE");
xlim([t0 tend])
grid on;
subplot(313);
plot(t, relposNEU(:,3))
title("relposU");
xlim([t0 tend])
grid on;

% figure(9);
% plot3(relvelNEU(:,2),relvelNEU(:,1),relvelNEU(:,3), '.');
% title("RELVEL")
% grid on;

figure(9);
plot3(relvelNEU(:,2),relvelNEU(:,1),relvelNEU(:,3), '.');
title("RELVEL")
grid on;

figure(99)
subplot(211)
plot(t, headMot), grid on
subplot(212)
plot(t, gSpeed), grid on, hold on
gSpeed_test = sqrt( relvelNEU(:,1).^2 + relvelNEU(:,2).^2 );
plot(t, gSpeed_test), hold off

figure(10);
subplot(311);
plot(t, relvelNEU(:,1))
title("relvelN");
xlim([t0 tend])
grid on;
subplot(312);
plot(t, relvelNEU(:,2))
title("relvelE");
xlim([t0 tend])
grid on;
subplot(313);
plot(t, relvelNEU(:,3))
title("relvelU");
xlim([t0 tend])
grid on;


%%
figure(234);
plot(t,DOP);
title("DOP")
legend(["gDOP" "pDOP" "tDOP" "vDOP" "hDOP" "nDOP" "eDOP"], "location", "northwest")
xlim([t0 tend])
grid on;

figure(124);
plot(t, K_pos);
title("Covariance Position")
grid on
legend(["NN" "NE" "ND" "EE" "ED" "DD"], 'location', 'best');
xlim([t0 tend])

figure(46);
plot(t,relaccNEU);
title("RELacc")
xlim([t0 tend])
grid on;

%%
figure(84)
plot(hpllh(:,1), hpllh(:,2), '.');
hold on
plot(llh(:,1), llh(:,2), '.');
hold off
title("hpposllh")
grid on;

figure(420)
% polarplot(heading, gSpeed,'.')
plot(t, headMot)
title("Ground Speed")
grid on;
%%
figure(456);
plot(t,D-[0 2 4 6 8 10 12 14])
ylim([-14 1])
legend(["invaldLLH" "RTK Fix" "RTK float" "diff available" "velCov valid" "posCov valid" "svin valid" "time mode"], "location", "southeast")
title("status")
xlim([t0 tend])
grid on;


