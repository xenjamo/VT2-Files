clc, clear variables
addpath ..\99_fcn_bib\
%%

% // rMat = Rz(-yaw) * Ry(pitch) * Rx (roll)
% // usual rpy euler angels use the convention Rz(yaw) * Ry(pitch) * Rx (roll) with the z axis pointing up (building a right hand system)
% int16_t roll;  /** roll lies in (-180, 180) deg
%                 *  rotating it 360 deg around x looks like this: 360      /   roll
%                 *                                                        /     ^
%                 *                                                       /      |
%                 *                                                  0   /        --> actual rotation */
% 
% int16_t pitch; /** pitch lies in (-90, 90) deg
%                 *  rotating it 360 deg around y looks like this:              pitch
%                 *                                                90    /\      ^     - this signal is dangerous for control
%                 *                                               -90      \/    |     - rpy euler angels have a singularity at -90 and +90 deg pitch
%                 *                                                               --> actual rotation */
% 
% int16_t yaw;   /** yaw lies in (0, 360) deg
%                 *  rotating it 360 deg around z looks like this: 3600      /   yaw
%                 *                                                         /     ^
%                 *                                                        /      |
%                 *                                                   0   /        --> actual rotation */
% 
% // DCM either transforms a vector from body to earth frame, e.g. Ev        = rMat * Bv  ( <->  Bv        = rMat^T * Ev )
% //         or rotates    a vector in the earth frame,       e.g. EvRotated = rMat * Ev  ( <->  BvRotated = rMat^T * Bv )
% // DCM  contains in the coloums [EexB  , EeyB  , EezB  ]  ( basis of the body  frame w.r.t. the earth frame )
% //      and      in the rows    [BexE^T; BeyE^T; BezE^T]  ( basis of the earth frame w.r.t. the body  frame )
% //  - without a mag the mahony filter yaw estimate has drift
% //  - the third row of the DCM (BezE^T) is invariant to yaw, this is why the mahony filter can estimate roll and pitch
% //    by only using an acc and a gyro
% //  - the static angle offset of roll and pitch which may bee seen after rearming are changes in accx and y due to temerature drift. 
% //    without an additional inertial sensor like a gnss unit it is not possible to estimate this drift

%%

% data = readmatrix('putty_14.log'); % gyro and acc
data = readmatrix('putty_15.log'); % gyro, acc and mag

time = data(:,10) * 1e-3;
time = time - time(1);
Ts = median(diff(time));
ind_gyro = 1:3;
ind_acc  = 4:6;
ind_mag  = 7:9;
ind_quat = 11:14;
ind_rpy  = 15:17;

figure(1)
plot(diff(time)), grid on

figure(2)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(ax(1), time, data(:,ind_gyro) * 180/pi), grid on, ylabel('Gyro (deg/sec)')
ax(2) = nexttile;
plot(ax(2), time, data(:,ind_acc)), grid on, ylabel('Acc (m/s^2)')
ax(3) = nexttile;
plot(ax(3), time, data(:,ind_mag)), grid on, ylabel('Mag'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

figure(3)
plot(time, cumtrapz(time, data(:,ind_gyro)) * 180/pi), grid on
ylabel('Gyro Integral (deg)'), xlabel('Time (sec)')
xlim([0, max(time)])

% bessel
p = 2;         % pole at p rad/s
kp = 2 * [p, p, p].';
ki = kp.^2 / 3;

% p = 7;         % pole at p rad/s
% kp = 2 * [p, p, p].';
% ki = kp.^2 / 3 * 0;

para.kp = kp;
para.ki = ki;

rpy0 = 0 * [60, 60, 60] * pi/180;
quat0 = rpy2quat(rpy0).';

[quatRP , biasRP ] = mahonyRP (data(:,ind_gyro), data(:,ind_acc), para, Ts, quat0);
[quatRPY, biasRPY] = mahonyRPY(data(:,ind_gyro), data(:,ind_acc), data(:,ind_mag), para, Ts, quat0);

rpyRP  = quat2rpy(quatRP );
rpyRPY = quat2rpy(quatRPY);

angleFun = @wrapToPi;

figure(4)
tiledlayout_ = tiledlayout(3,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
ang_ind = 1;
plot(time, angleFun(cumtrapz(time, data(:,ind_gyro(ang_ind)))) * 180/pi), grid on, hold on
plot(time, angleFun(rpyRP (:,ang_ind)) * 180/pi)
plot(time, angleFun(rpyRPY(:,ang_ind)) * 180/pi)
plot(time, angleFun(data(:,ind_rpy(ang_ind))) * 180/pi), hold off
ylabel('Roll (deg)'), xlabel('Time (sec)')
ax(2) = nexttile;
ang_ind = 2;
plot(time, angleFun(cumtrapz(time, data(:,ind_gyro(ang_ind)))) * 180/pi), grid on, hold on
plot(time, angleFun(rpyRP (:,ang_ind)) * 180/pi)
plot(time, angleFun(rpyRPY(:,ang_ind)) * 180/pi)
plot(time, angleFun(data(:,ind_rpy(ang_ind))) * 180/pi), hold off
ylabel('Pitch (deg)'), xlabel('Time (sec)')
ax(3) = nexttile;
ang_ind = 3;
plot(time, angleFun(cumtrapz(time, data(:,ind_gyro(ang_ind)))) * 180/pi), grid on, hold on
plot(time, angleFun(rpyRP (:,ang_ind)) * 180/pi)
plot(time, angleFun(rpyRPY(:,ang_ind)) * 180/pi)
plot(time, angleFun(data(:,ind_rpy(ang_ind))) * 180/pi), hold off
ylabel('Yaw (deg)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

figure(5)
tiledlayout_ = tiledlayout(2,1); tiledlayout_.TileSpacing = 'compact';
ax(1) = nexttile;
plot(time, biasRP  * 180/pi), grid on
ylabel('Gyro Bias (deg/s)')
ax(2) = nexttile;
plot(time, biasRPY * 180/pi), grid on
ylabel('Gyro Bias (deg/s)'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax, xlim([0, max(time)])

% figure(6)
% subplot(211)
% plot(time, data(:,ind_quat)), grid on
% ylabel('Quaternion')
% subplot(212)
% plot(time, quat_), grid on
% ylabel('Quaternion'), xlabel('Time (sec)')
