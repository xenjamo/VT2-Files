fid = fopen('acc_calib.bin');
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
    M = [M;M_];                     % write array to matrix
    A = [A;A_];
    A2 = [A2;A2_];
    B = [B;B_];
    C = [C;C_];
    D = [D;bitget(D_,8:-1:1)];
end

fclose(fid);                        % close file

%%
close all;
size(M)
B(end,2)
t = B(:,1)/1000000;
dt = [0;diff(t)];
t = t-t(1);
figure(999);
plot(t,dt);
title("perf")
grid on;

figure(1);
plot(t,A(:,1));
hold on;
plot(t,A(:,2));
plot(t,A(:,3));
hold off;
title("gyro")

figure(2);
hold on;
plot(t,A(:,6));
hold off;
title("acc")

figure(3);
plot(t,A(:,7));
hold on;
plot(t,A(:,8));
plot(t,A(:,9));
hold off;
title("mag")

figure(5);
plot(t,A(:,14));
hold on;
plot(t,A(:,15));
plot(t,A(:,16));
hold off;
title("rpy")


figure(234);
hold on;
plot(t,A(:,30));
plot(t,A(:,31));
plot(t,A(:,32));
plot(t,A(:,33));
plot(t,A(:,34));
plot(t,A(:,35));
plot(t,A(:,36));
hold off;
title("DOP")

figure(2674);
hold on;
plot3(A(:,37),A(:,38),-A(:,39), '.');
hold off;
title("RELPOS")

figure(46);
hold on;
plot(t,A(:,40));
plot(t,A(:,41));
plot(t,A(:,42));
hold off;
title("RELacc")

figure(84)
hold on;
plot3(A2(:,1), A2(:,2), A2(:,3), '.');
hold off;
title("hpposllh")

figure(456);
plot(t,D-[0 2 4 6 8 10 12 14])
ylim([-14 1])
legend({"invaldLLH" "RTK Fix" "RTK float" "diff available" "velCov valid" "posCov valid" "svin valid" "time mode"}, "location", "southeast")
title("status")
