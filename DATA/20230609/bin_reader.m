clear all, clc, close all
%% INDEX
loc = 'bin_files/';
% filename = '001.bin' % under bridge;
% filename = '002.bin' % under bridge;
% filename = '003.bin' % under bridge;
% filename = '004.bin' % broken;
% filename = '005.bin' % mag calib file;
% filename = '006.bin' % general flight;
% filename = '007.bin' % line walk;
% filename = '008.bin' % line walk;
% filename = '011.bin' % Roll Pitch Yaw;
% filename = '013.bin' % speed;
% filename = '014.bin' % speeeed;
% filename = '015.bin' % XYZ;
% filename = '016.bin' % XYZ;
% filename = '017.bin' % chirp;
% filename = '018.bin' % acro;
filename = '020.bin' % acro;
% filename = '021.bin' % acro;
% filename = '022.bin' % acro;




%% read file

[A,A2,B,C,D] = decode_bin(append(loc,filename));
M = [A,A2,B,C,D];

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

A_mag = [ 0.9821680,  0.0000000,  0.0000000;
          0.0159774,  0.9866579,  0.0000000;
          0.0536142,  0.0066898,  1.0311741];
b_mag = [-0.4711387,  0.3239450, -0.3728509];
mag_calib = (mag - b_mag)*A_mag.';

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

hpllh = A2(:,1:3);

save data.mat n t itow ttff msss rtk_flags...
     gnss_fix fix_type numSV lastCorrectionAge...
     gyro accel mag quat rpy K_pos K_vel...
     DOP relposNED relaccNED hMSL hAcc vAcc...
     relvelNED sAcc gSpeed headMot headAcc...
     magDec magAcc hpllh mag_calib 

%% plots

figure(999);
subplot(411)
plot(t,gnss_fix);
title("GPS Fix")
grid on;
subplot(412)
plot(t,rtk_flags);
title("rtk flags")
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

figure(87);
subplot(211)
plot(t, sqrt(sum(mag.^2, 2)));
title("mag")
xlim([t0, tend])
grid on;
subplot(212);
plot(t, sqrt(sum(mag_calib.^2, 2)));
title("mag calib")
xlim([t0, tend])
grid on;

% carefull here
% Tmean = 6.0;
% ind0 = t > 1 & t < Tmean;
% mag0 = median(mag_calib(ind0)).';
% Ymag0 = atan2(mag_calib(ind0, 2), mag_calib(ind0, 1));
% Ymag0 = median(Ymag0);
Ymag = atan2(mag(:,2), mag(:,1));
% Ymag = atan2(sin(Ymag - Ymag0), cos(Ymag - Ymag0)); % unwrap(Ymag - Ymag0);
figure(33)
plot(t, [Ymag, headMot]), grid on

figure(5);
plot(t,rpy);
title("rpy")
xlim([t0, tend])
grid on;
%%
figure(7);
plot(relposNED(:,2),relposNED(:,1), '.-');
title("RELPOS")
grid on;

figure(8);
subplot(311);
plot(t, relposNED(:,1))
title("relposN");
xlim([t0 tend])
grid on;
subplot(312);
plot(t, relposNED(:,2))
title("relposE");
xlim([t0 tend])
grid on;
subplot(313);
plot(t, relposNED(:,3))
title("relposD");
xlim([t0 tend])
grid on;

% figure(9);
% plot3(relvelNEU(:,2),relvelNEU(:,1),relvelNEU(:,3), '.');
% title("RELVEL")
% grid on;
%% 
figure(9);
plot3(relvelNED(:,2),relvelNED(:,1),relvelNED(:,3), '.');
title("RELVEL")
grid on;

figure(90)
subplot(211)
plot(t, headMot); grid on
title("headMot")
subplot(212)
plot(t, gSpeed); grid on
title("gSpeed")

figure(901)
subplot(211)
plot(t, headAcc); grid on
title("headAcc")
subplot(212)
plot(t, gSpeed); grid on
title("gSpeed")


figure(10);
subplot(311);
plot(t, relvelNED(:,1))
title("relvelN");
xlim([t0 tend])
grid on;
subplot(312);
plot(t, relvelNED(:,2))
title("relvelE");
xlim([t0 tend])
grid on;
subplot(313);
plot(t, relvelNED(:,3))
title("relvelU");
xlim([t0 tend])
grid on;


%%
figure(11);
plot(t,DOP);
title("DOP")
legend(["gDOP" "pDOP" "tDOP" "vDOP" "hDOP" "nDOP" "eDOP"], "location", "northwest")
xlim([t0 tend])
grid on;

figure(12);
plot(t, K_pos);
title("Covariance Position")
grid on
legend(["NN" "NE" "ND" "EE" "ED" "DD"], 'location', 'best');
xlim([t0 tend])

figure(13);
plot(t,hAcc);
title("hAcc")
xlim([t0 tend])
grid on;

%%
figure(14)
plot(hpllh(:,1), hpllh(:,2), '.');
title("hpposllh")
grid on;

figure(140);
subplot(311);
plot(t,hpllh(:,1)); grid on;
ylabel('lon');
subplot(312);
plot(t,hpllh(:,2)); grid on;
ylabel('lat');
subplot(313)
plot(t,hAcc);


figure(15)
% polarplot(heading, gSpeed,'.')
plot(t, headMot)
title("headMot")
grid on;
%%
figure(16);
plot(t,D-[0 2 4 6 8 10 12 14])
ylim([-14 1])
legend(["invaldLLH" "RTK Fix" "RTK float" "diff available" "velCov valid" "posCov valid" "svin valid" "time mode"], "location", "southeast")
title("status")
xlim([t0 tend])
grid on;


