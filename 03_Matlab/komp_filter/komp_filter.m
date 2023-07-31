load data_rp.mat
%%
w = 3.14;
Gtp = tf([w], [1 w]);
Ghp = 1-Gtp;
tend = t(end);

%%
[m,n] = size(accel);

a_yz = [0 0 0];
a_xz = [0 0 0];


rpy = zeros(m,3);

%%
for i = 1:m
    an = accel(i,:) / norm(accel(i,:));
    a_yz(2:3) = an(2:3);
    a_xz([1,3]) = an([1,3]);

    [roll,vr] = calcAngleBetween2Vectors(a_yz', [0 0 1]');
    [pitch,vp] = calcAngleBetween2Vectors(a_xz', [0 0 1]');
    rpy(i,:) = [roll*norm(a_yz)*max(vr), pitch*norm(a_xz)*max(vp), 0];
end


figure(1)
plot(t,rpy);grid on;
title('RPY acc');

figure(3)
plot3(accel(:,1),accel(:,2),accel(:,3));grid on; axis equal;
xlim([-12 12]), ylim([-12 12]), zlim([-12 12]); 
title('RPY acc');


%%
figure(2)
bode(Gtp,Ghp);