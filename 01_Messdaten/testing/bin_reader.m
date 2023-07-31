fid = fopen('024.bin');
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
t = B(:,1)/1000;
dt = diff(t);
figure(666)
title("performance in dt [s]")
plot(dt)

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

figure(20);
title("bleb")
plot(t,A(:,42));
hold on;
plot(t,A(:,43));
plot(t,A(:,44));
hold off;



