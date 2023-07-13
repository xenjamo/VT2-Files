function [quat_, bias_] = mahonyRPY(gyro, acc, mag, para, Ts, quat0)

if nargin == 5
    quat = eye(4,1);
else
    quat = quat0;
end

kp = para.kp;
ki = para.ki;

N = size(gyro, 1);
bias = zeros(3,1);
quat_ = zeros(N, 4);
bias_ = zeros(N, 3);
for i = 1:N

    % h = (hx; hy;  0) - measured mag field vector in EF (assuming Z-component is zero)
    % b = (bx;  0;  0) - reference mag field vector heading due North in EF (assuming Z-component is zero)
    mag_n = mag(i,:).';
    mag_n = mag_n ./ norm(mag_n);

    CBE = quat2CBE(quat);
    CEB = CBE.';
    h = CEB * mag_n; % Cmu.' * CEB = (CEB.' * Cmu).' = = (CBE * Cmu).'
    h(3) = 0;

    % b = [1; 0; 0];
    % e = CBE * cross(h/norm(h), b);
    
    % by normalising here you adjust implicit the weigth
    b = [norm(h); 0; 0];
    e = CBE * cross(h, b); % e_heading = CBE * cross(h, b);

    g_n = CEB(3,:).';

    acc_n = acc(i,:).';
    acc_n = acc_n ./ norm(acc_n);
    e = e + cross(acc_n, g_n); % e_level = cross(acc_n, g_n)
        
    bias = bias + ki .* e * Ts;

    Q = [[-quat(2), -quat(3), -quat(4)]; ...
         [ quat(1), -quat(4),  quat(3)]; ...
         [ quat(4),  quat(1), -quat(2)]; ...
         [-quat(3),  quat(2),  quat(1)]];

    % dquat = Ts * 0.5 * Q * ( gyro(i,:).' + bias + kp_level .* e_level + kp_heading .* e_heading);
    dquat = Ts * 0.5 * Q * ( gyro(i,:).' + bias + kp .* e );
    quat = quat + dquat;
    quat = quat ./ norm(quat);

    quat_(i,:) = quat.';
    bias_(i,:) = bias.';
end


return