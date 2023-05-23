fid = fopen('001.bin'); % this file format is [3x Gyro_float, 3x Acc_float, 3x Mag_float]


A = fread(fid,'float');


fclose(fid);

% [fm, fn] = size(A);
% m = 11;
% n = fm/m;
% 
% A = reshape(A,[m,n])'

