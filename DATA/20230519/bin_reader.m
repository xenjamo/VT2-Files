fid = fopen('001.bin');
M = []; % initalize array
A = [];
B = [];
C = [];
D = [];

while 1
    A_ = fread(fid,[1 50],'float');
    B_ = fread(fid,[1 6],'uint32');
    C_ = fread(fid,[1 4],'uint8');
    D_ = fread(fid,[1 8],'uint8');
    M_ = [A_,B_,C_,D_];                   
    if size(M_) < 68                
        break;                      
    end
    M = [M;M_];                     % write array to matrix
    A = [A;A_];
    B = [B;B_];
    C = [C;C_];
    D = [D;D_];
end

fclose(fid);                        % close file
%% 
size(M)
B(end,2)
t = B(:,1)/1000;
dt = [0;diff(t)];
figure(999);
title("perf")
plot(t,dt);
grid on;

figure(1);
title("gyro")
plot(t,A(:,1));
hold on;
plot(t,A(:,2));
plot(t,A(:,3));
hold off;

figure(2);
title("acc")
plot(t,A(:,4));
hold on;
plot(t,A(:,5));
plot(t,A(:,6));
hold off;

figure(234);
title("DOP")
hold on;
plot(t,A(:,23));
plot(t,A(:,24));
plot(t,A(:,25));
plot(t,A(:,26));
plot(t,A(:,27));
plot(t,A(:,28));
plot(t,A(:,29));
hold off;

figure(2674);
title("RELPOS")
hold on;
plot(t,A(:,30));
plot(t,A(:,31));
plot(t,A(:,32));
hold off;

figure(456);
title("RELacc")
hold on;
plot(t,A(:,33));
plot(t,A(:,34));
plot(t,A(:,35));
hold off;

figure(20);
title("acc")
plot(t,A(:,42));
hold on;
plot(t,A(:,43));
plot(t,A(:,44));
hold off;


figure(9865);
title("status");
subplot(311);
plot(t,D(:,5))
ylabel("RTCM")
subplot(312);
plot(t,D(:,6))
ylabel("RTK float")
subplot(313);
plot(t,D(:,7))
ylabel("RTK fix")



