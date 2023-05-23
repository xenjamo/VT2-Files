fid = fopen('011.bin');

A = fread(fid,[1 11],'fffffffff ii');



fclose(fid);

% [fm, fn] = size(A);
% m = 11;
% n = fm/m;
% 
% A = reshape(A,[m,n])'

