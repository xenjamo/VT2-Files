function [A,A2,B,C,D] = decode_bin(filename)
%Decodes a binary file from the VT2 Project

fid = fopen(filename);
M = []; % initalize array
A = [];
A2 = [];
B = [];
C = [];
D = [];

while 1
    A_ = fread(fid,[1 54],'float');
    A2_ = fread(fid,[1 3], 'double');
    B_ = fread(fid,[1 6],'uint32');
    C_ = fread(fid,[1 4],'uint8');
    D_ = fread(fid,[1 1],'uint8');
    M_ = [A_,A2_,B_,C_,D_];                   
    if size(M_) < 68                
        break;                      
    end
    % M = [M;M_];                     % write array to matrix
    A = [A;A_];
    A2 = [A2;A2_];
    B = [B;B_];
    C = [C;C_];
    D = [D;bitget(D_,8:-1:1)];
end
M = [A,A2,B,C,D];

fclose(fid);                        % close file

end

