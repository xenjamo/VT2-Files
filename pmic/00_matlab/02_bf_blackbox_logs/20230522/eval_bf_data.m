clc, clear variables
addpath ..\..\99_fcn_bib\
%%

% was_chirp_excitation = true;
% 
% fileName = 'LOG00002.BFL.csv';
% [para, NHeader] = get_para_and_NHeader(fileName);
% 
% tic
% % read the data
% data = readmatrix(fileName, 'NumHeaderLines', NHeader);
% % data = data(2*2e3:end,:); % in case of logging issues and you want to neglect some data at the beginning
% [~, Nsig] = size(data)
% 
% 
% % check if there are NaN entries % ToDo: why are there nans from col 49-56
% ind_nan = max(sum(isnan(data)));
% data = data(ind_nan+1:end,:);
% [Ndata, Nsig] = size(data);
% if max(sum(isnan(data))) ~= 0
%     error('There are still NaN entries in the data matrix.')
% end
% toc
% 
% % --- HACK STARTS ---
% % data = data(178*2e3:end,:);
% % Ndata = size(data, 1);
% % --- HACK ENDS   ---
% 
% % eval time
% time = (data(:,2) - data(1,2))*1e-6;
% dtime_meas_mus = diff(time)*1e6;
% 
% figure(1)
% stairs(time(1:end-1), dtime_meas_mus), grid on
% title(sprintf('mean: %0.3f mus, median: %0.3f mus, std: %0.3f mus\n', ...
%     mean(dtime_meas_mus), ...
%     median(dtime_meas_mus), ...
%     std(dtime_meas_mus)))
% xlabel('Time (sec)'), ylabel('Ts log (mus)')
% xlim([0, max(time)])

% create signal indexes % ToDo: change this logic
ind_pid_p   =  3: 5;
ind_pid_i   =  6: 8;
ind_pid_d   =  9:10;
ind_pid_f   = 11:13;
ind_rc      = 14:17;
ind_setp    = 18:21;
ind_baro    = 24;
ind_gyro_f  = 26:28;
ind_acc     = 29:31;
ind_mot     = 40:43;
ind_gyro    = 32:34;
ind_motors  = 36:39;
ind_sinarg  = 35;
ind_pid_sum = 58:60;
ind_pid_err = 65:67;
ind_pi_sum  = 72:74; % will be created below

% % unscaling from highResolutionGain
% if para.blackbox_high_resolution
%     blackboxHighResolutionScale = 10.0;
%     ind_bb_high_res = [ind_gyro_f, ind_rc, ind_setp(1:3)];
%     data(:, ind_bb_high_res) = 1/blackboxHighResolutionScale * data(:, ind_bb_high_res);
% end
% 
% % unscaling and remapping
% if was_chirp_excitation
%     gyroScale = 1.6e1;
%     data(:, ind_gyro) = 1/gyroScale * data(:, ind_gyro);
%     sinargScale = 5.0e3;
%     data(:, ind_sinarg) = 1/sinargScale * data(:, ind_sinarg);
% end
% 
% % % reconstructing, original data will be overwritten
% % data(:,ind_pid_sum) = data(:,ind_pid_p) + data(:,ind_pid_i) + [data(:,ind_pid_d), zeros(Ndata,1)] + data(:,ind_pid_f);
% % data(:,ind_pid_err) = data(:,ind_setp(1:3)) - data(:,ind_gyro_f);
% 
% % negative sign for pid_err
% data(:,ind_pid_err) = -data(:,ind_pid_err);
% 
% % create an additional entry for the pi sum
% data = [data, data(:,ind_pid_p) + data(:,ind_pid_i)];
% 
% % different rates
% Ts      = para.looptime * 1e-6;             % gyro (base loop)
% Ts_cntr = para.pid_process_denom*Ts;        % cntrl
% Ts_log  = para.frameIntervalPDenom*Ts_cntr; % logging

load LOG00002_bfl.mat % save LOG00002_bfl data Ts Ts_cntr Ts_log time

%% show gyro to select Teval and spectras (gyro and pid sum)

figure(2)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), time, data(:, [ind_setp(1), ind_gyro(1), ind_gyro_f(1)])), grid on, ylabel('Roll (deg/sec)')
ax(2) = nexttile;
plot(ax(2), time, data(:, [ind_setp(2), ind_gyro(2), ind_gyro_f(2)])), grid on, ylabel('Pitch (deg/sec)')
ax(3) = nexttile;
plot(ax(3), time, data(:, [ind_setp(3), ind_gyro(3), ind_gyro_f(3)])), grid on, ylabel('Yaw (deg/sec)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

data_spectras = data(:, [ind_gyro, ind_gyro_f, ind_pid_sum, ind_setp(1:3)]);

% Nest     = round(2.5/Ts_log);
% koverlap = 0.95;
% Noverlap = round(koverlap*Nest);
% window   = hann(Nest);
% [pxx, f] = pwelch(data_spectras, window, Noverlap, Nest, 1/Ts_log, 'power');

[pxx, f] = pwelch(data_spectras, round(10/Ts_log), [], [], 1/Ts_log, 'Power');

spectras = sqrt(pxx*2); % power -> amplitude (dc needs to be scaled differently)

figure(3)
tiledlayout_ = tiledlayout(2,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
ind_show = 1:6;
plot(ax(1), f, spectras(:, ind_show)), grid on, set(gca, 'YScale', 'log'), ylabel('Gyros (deg/sec)')
ax(2) = nexttile;
ind_show = 7:9;
plot(ax(2), f, spectras(:, ind_show)), grid on, set(gca, 'YScale', 'log'), ylabel('PID Sum')
% ax(3) = nexttile;
% ind_show = 10:12;
% plot(ax(3), f, spectras(:, ind_show)), grid on, set(gca, 'YScale', 'log'), ylabel('Setpoint (deg/sec)'), xlabel('Frequency (Hz)')
linkaxes(ax), clear ax, axis([0 1/2/Ts_log 1e-3 1e1])

figure(4)
tiledlayout_ = tiledlayout(2,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), time, data(:, ind_setp(1:3))), grid on, ylabel('Setpoint (deg/sec)')
ax(2) = nexttile;
plot(ax(2), time, data(:, ind_sinarg)), grid on, ylabel('Sinarg (rad)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

figure(5)
tiledlayout_ = tiledlayout(4,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), time, data(:, ind_gyro)), grid on, ylabel('GyroS (deg/sec)')
ax(2) = nexttile;
plot(ax(2), time, data(:, ind_pid_sum)), grid on, ylabel('PID Sum')
ax(3) = nexttile;
plot(ax(3), time, data(:, ind_mot)), grid on, ylabel('Motors')
ax(4) = nexttile;
plot(ax(4), time, data(:, ind_setp(4))), grid on, ylabel('Throttle'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax
xlim([0, max(time)])

figure(6)
plot(time, data(:,ind_acc)), grid on

figure(7)
plot(time, data(:,ind_gyro)), grid on


