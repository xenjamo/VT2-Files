clc, clear variables
addpath ..\..\99_fcn_bib
%%

% gyro mag time
load data_20230609_pb_mag_020.mat
% load data_20230609_pb_mag_021.mat
% load data_20230609_pb_mag_022.mat

% apply calibration
load 20230609_pb_mag_calib_005.mat % a_mag b_mag
mag  = (mag - b_mag.') * A_mag.';

mat_name = 'data_20230609_bf_020';
% mat_name = 'data_20230609_bf_021';
% mat_name = 'data_20230609_bf_022';

[ind, gyro_bf, acc_bf, mag_bf] = get_ind_and_synched_bf_data_20230609(gyro, mat_name);
gyro_pb = gyro(ind, :);
mag_pb  = mag(ind, :);

figure(1)
plot([gyro_pb,  gyro_bf])

figure(2)
subplot(211)
plot( [mag_pb, mag_bf] )
subplot(212)
plot( [sqrt( sum( mag_pb.^2, 2) ), sqrt( sum( mag_bf.^2, 2) )] )
