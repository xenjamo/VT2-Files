function [ind, gyro_bf, acc_bf, mag_bf] = get_ind_and_synched_bf_data_20230609(gyro, mat_name)
% it is assumed that pes board was always running before and after bf bbl

% parameters
mag_alignment_roll  = 180;
mag_alignment_pitch = -29;
fg = 10;

% load mat
eval(['load ', mat_name, '.mat']);
ind_gyro = ind_gyro_f; clear ind_gyro_f

% bf data preprocessing
data(:, ind_gyro) = data(:, ind_gyro) * pi / 180;    % rad/sec

data(:, ind_acc ) = data(:, ind_acc ) / 2000 * 9.81; % m/s^2

data(:, ind_mag)  = data(:, ind_mag) / 2000;
R = calcRn([0 1 0].', mag_alignment_pitch * pi/180) * ...
    calcRn([1 0 0].', mag_alignment_roll  * pi/180);
data(:, ind_mag) = data(:, ind_mag) * R.';


% apply bf calibration
load 20230609_bf_mag_calib_mean.mat
data(:, ind_mag)  = (data(:, ind_mag) - b_mag_bf.') * A_mag_bf.';


% filtering
Gf = c2d( tf(1, [1/(2*pi*fg) 1] ), Ts_log, 'tustin');
data(:, [ind_gyro, ind_acc, ind_mag]) = ...
    filtfilt(Gf.num{1}, Gf.den{1}, data(:, [ind_gyro, ind_acc, ind_mag]));


% downsampling
n_ds = (1/50) / (1/2000);
data = data(1:n_ds:end, :);

% cut gyro to to bf length
N_bf = size(data(:,1), 1);
gyro = gyro(1:N_bf, :);

Gf = c2d( tf(1, [1/(2*pi*fg) 1] ), 1/50, 'tustin');
gyro_corr = xcorr(filtfilt(Gf.num{1}, Gf.den{1}, gyro(:,1)), ...
                  data(:, ind_gyro(1)), ...
                  'coeff');

% shift bf data accordingly
[~, ind_corr_max] = max(gyro_corr);
n_shift = N_bf - ind_corr_max;
ind = false( size(gyro, 1), 1);
if n_shift > 0
    ind(1:end-n_shift, :) = true;
    data = data(n_shift+1:end, :);
end

gyro_bf = data(:, ind_gyro);
acc_bf  = data(:, ind_acc );
mag_bf  = data(:, ind_mag );

return
