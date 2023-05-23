%%

clc, clear variables
%%
fid = fopen('acc_calib.bin');
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

acc = data(1:350,4:6);
acc_mean = mean(acc);
b_acc = [acc_mean(:,1:2),0];
acc_calib = acc - b_acc;

figure(1)
ax(1) = subplot(221);
plot(acc), grid on
title('uncalibrated')
ax(2) = subplot(223);
plot(sqrt(sum(acc.^2, 2))), grid on
ax(3) = subplot(222);
plot(acc_calib), grid on
title('calibrated')
ax(4) = subplot(224);
plot(sqrt(sum(acc_calib.^2, 2))), grid on
linkaxes(ax, 'x'), clear ax
xlim([0 size(acc, 1)])

fprintf('b_acc << %10.7ff, %10.7ff, %10.7ff;\n', b_acc)
