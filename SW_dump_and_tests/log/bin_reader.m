fid = fopen('011.bin');
M = []; % initalize array

while 1
    A = fread(fid,[1 9],'float');   % reads 9x a float value and puts them in a 9x1 array
    B = fread(fid,[1 2],'int');     % reads 2x a int value and puts them in a 2x1 array
    M_ = [A,B];                     % merge both to a single array (not necessary can be kept separatly)
    if size(M_) == 0                % check whether this array is empty
        break;                      % if its empty we have reached EOF
    end
    M = [M;M_];                     % write array to matrix
end

fclose(fid);                        % close file
