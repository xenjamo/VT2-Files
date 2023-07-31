function [t] = saveForLatex1d(var,time,filename,margin)
t = 1;

ymid = (max(max(var)) + min(min(var))) / 2;
dy = (max(max(var)) - min(min(var))) / 2;

y_lim = [-dy*margin, dy*margin] + ymid;
writematrix([time,var], append('csv/',filename));
filename
fprintf('xmin=%.7f, xmax=%.7f,\nymin=%.7f, ymax=%.7f,\n', time(1), time(end),...
                                                           y_lim(1), y_lim(2));

end