fid = fopen('015.bin');
M = []; % initalize array
A = [];
B = [];
C = [];
D = [];

while 1
    A_ = fread(fid,[1 50],'float');
    B_ = fread(fid,[1 6],'uint32');
    C_ = fread(fid,[1 4],'uint8');
    D_ = fread(fid,[1 8],'uint8');
    M_ = [A_,B_,C_,D_];                   
    if size(M_) < 68                
        break;                      
    end
    M = [M;M_];                     % write array to matrix
    A = [A;A_];
    B = [B;B_];
    C = [C;C_];
    D = [D;D_];
end

fclose(fid);                        % close file
%% 
size(M)


figure(1);
plot(A(:,14));

