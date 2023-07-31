clc, clear variables
%%

% https://github.zhaw.ch/altb/IndNav_RPI_computer.git
% SHA-1: 0f7c8f1

do_export_as_pdfs = false;

%% session0.log

mission_number = 1; % 1, 2, 4

% \src\Mapper.cpp
%     m_data->Printf("%e %e %e %e %e %e\n", m_pos_i.x(), m_pos_i.y(), m_pos_i.z(), rpy_i.x(), rpy_i.y(), rpy_i.z());

file_path = ['..\00_Data\2nd_floor_RPI4_data\log_LS_3rd_floor_', ... 
             num2str(mission_number), ...
             '\DATA\session0.log'];
data_table = readtable( file_path );
data_mapper = table2array( data_table(:, 2:end) );

% time resolution sec is not usefull
N_data = size(data_mapper, 1);
time_mapper = zeros(N_data, 1);
for i = 1:N_data
    time_str = data_table{i,1}{1}(end-5:end);
    time_mapper(i) = str2double( time_str(1:2) ) * 60 * 60 + ...
                     str2double( time_str(3:4) ) * 60 + ...
                     str2double( time_str(5:6) );
end
time_mapper = time_mapper - time_mapper(1);

% \src\FlirCntrl.cpp: Line 160 - 162
%     m_logger_cam_pose->Printf("%e %e %e %e %e %e\n",
%                               MAP_est_xyz(0), MAP_est_xyz(1), MAP_est_xyz(2),
%                               MAP_est_RPY(0), MAP_est_RPY(1), MAP_est_RPY(2));
file_path = ['..\00_Data\2nd_floor_RPI4_data\log_LS_3rd_floor_', ... 
             num2str(mission_number), ...
             '\FLIR_POSE\session0.log'];
data_table = readtable( file_path );
data_cam_pose = table2array( data_table(:, 2:end) );

lw = 2;

pdfsize = [];

ind_pos = 1:3;
ind_rpy = 4:6;

rgb_colors = [[1 0   0]; ...
              [0 0.5 0]; ...
              [0 0   1]];

figure(1)
plot( time_mapper, data_mapper(:, ind_pos)), grid on
xlabel('Time (sec)'), ylabel('Position (m)')
legend('x', 'y', 'z', 'location', 'best')
set(findall(gcf, 'type', 'line'), 'linewidth', lw)

figure(2)
plot( time_mapper, data_mapper(:, ind_rpy) * 180/pi), grid on
xlabel('Time (sec)'), ylabel('Angles (deg)'), legend
legend('roll', 'pitch', 'yaw', 'location', 'best')
set(findall(gcf, 'type', 'line'), 'linewidth', lw)

figure(3)
plot3( data_mapper(:, ind_pos(1)), data_mapper(:, ind_pos(2)), data_mapper(:, ind_pos(3)) ), hold on
xlabel('x-Axis (m)'), ylabel('y-Axis (m)'), zlabel('z-Axis (m)')
for i = 1:size(data_cam_pose, 1)
    pos_i = data_cam_pose(i, ind_pos);
    R_i = eul2rotm( fliplr( data_cam_pose(i, ind_rpy) ), 'ZYX');
    arrow_length = 0.2;
    R_i_scaled = R_i * arrow_length;
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,1), R_i_scaled(2,1), R_i_scaled(3,1), 'LineWidth', 2, 'AutoScale', 'off', 'color', [1 0 0])
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,2), R_i_scaled(2,2), R_i_scaled(3,2), 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0.5 0])
    quiver3(pos_i(1,1), pos_i(1,2), pos_i(1,3), R_i_scaled(1,3), R_i_scaled(2,3), R_i_scaled(3,3), 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0 1])
end
hold off, grid on, axis equal, zlim([0 0.6])
set(findall(gcf, 'type', 'line'), 'linewidth', lw)

if do_export_as_pdfs
    pdfexport(['flir_pose_mission_', num2str(mission_number)], ...
              [], ...
              pdfsize);
end

%% wp_log.txt

file_path = ['..\00_Data\2nd_floor_RPI4_data\log_LS_3rd_floor_', ... 
             num2str(mission_number), ...
             '\wp_log.txt'];
data_wp = readmatrix( file_path );

lw = 2;

pdfsize = 0.8*[15 12];

ind_pos_des = 5:7;
ind_pos = 2:4;
ind_yaw_des = 9;
ind_yaw = 8;

N_data = size(data_wp, 1);
time_wp = (0:N_data-1) * time_mapper(end) / (N_data - 1);

figure(4)
plot( time_wp, data_wp(:, [ind_pos_des, ind_pos]) ), grid on
xlabel('Time (sec)'), ylabel('Position (m)')
xlim([0 time_wp(end)])
legend('x', 'y', 'z', 'location', 'best')
set(findall(gcf, 'type', 'line'), 'linewidth', lw)

if do_export_as_pdfs
    pdfexport(['tracking_comparison_position_mission_', num2str(mission_number)], ...
              [], ...
              pdfsize);
end

figure(5)
plot( time_wp, data_wp(:, [ind_yaw_des, ind_yaw]) ), grid on
xlabel('Time (sec)'), ylabel('Yaw (deg)')
xlim([0 time_wp(end)])
set(findall(gcf, 'type', 'line'), 'linewidth', lw)

if do_export_as_pdfs
    pdfexport(['tracking_comparison_yaw_mission_', num2str(mission_number)], ...
              [], ...
              pdfsize);
end


