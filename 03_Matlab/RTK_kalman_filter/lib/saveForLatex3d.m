function [t] = saveForLatex3d(var,filename,margin)
t = 1;
xmid = (max(var(:,1)) + min(var(:,1))) / 2;
dx = (max(var(:,1)) - min(var(:,1))) / 2;
ymid = (max(var(:,2)) + min(var(:,2))) / 2;
dy = (max(var(:,2)) - min(var(:,2))) / 2;
zmid = (max(var(:,3)) + min(var(:,3))) / 2;
dz = (max(var(:,3)) - min(var(:,3))) / 2;

x_lim = [-dx*margin, dx*margin] + xmid;
y_lim = [-dy*margin, dy*margin] + ymid;
z_lim = [-dz*margin, dz*margin] + zmid;
filename
writematrix(var, filename);
fprintf('xmin=%.7f, xmax=%.7f,\nymin=%.7f, ymax=%.7f,\nzmin=%.7f, zmax=%.7f,\n',...
                                                           x_lim(1), x_lim(2),...
                                                           y_lim(1), y_lim(2),...
                                                           z_lim(1), z_lim(2));
end

