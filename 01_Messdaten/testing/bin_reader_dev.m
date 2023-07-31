%% make sure ur using the correct format


%% read file
fid = fopen('007_.bin');
M = []; % initalize array
A = [];
A2 = [];
B = [];
C = [];
D = [];

while 1
    A_ = fread(fid,[1 54],'float');
    A2_ = fread(fid,[1 3],'double'); % if no double format is read change size to [1 0]
    B_ = fread(fid,[1 6],'uint32');
    C_ = fread(fid,[1 4],'uint8');
    D_ = fread(fid,[1 1],'uint8'); %may increase
    M_ = [A_,A2_,B_,C_,D_];                   
    if size(M_) < 68                
        break;                      
    end
    M = [M;M_];                     % write array to matrix
    A = [A;A_];
    A2 = [A2;A2_];
    B = [B;B_];
    C = [C;C_];
    D = [D;D_];
end

fclose(fid);                        % close file
%% decode bits
[m,n] = size(D);
d = [];
for i = 1:m
    d_ = bitget(D(i),1:1:8);
    d = [d;d_];
end
%% 
size(M)
t = B(:,1)/1000000;
t = t - t(1);
dt = diff(t);
figure(666)
plot(dt)
title("performance in dt [s]")
ylim([0 0.03])
grid on;

figure(1);
plot(t,A(:,1));
hold on;
plot(t,A(:,2));
plot(t,A(:,3));
hold off;
title("gyro")

figure(2);
plot(t,A(:,4));
hold on;
plot(t,A(:,5));
plot(t,A(:,6));
hold off;
title("acc")

figure(5);
plot(t,A(:,14));
hold on;
plot(t,A(:,15));
plot(t,A(:,16));
hold off;
title("rpy")

if A2
    figure(20);
    plot(t,A2(:,1));
    hold on;
    plot(t,A2(:,2));
    plot(t,A2(:,3));
    hold off;
    title("bleb")
end

%% status
figure(4658)
dp = d(:,:) + [0 2 4 6 8 10 12 14];
plot(t,dp(:,:));
grid on;
set(gca,'ytick',[0:1:14])



