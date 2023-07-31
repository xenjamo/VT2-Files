function [para, NHeader] = get_para_and_NHeader(filePath)

para = [];

fid = fopen(filePath);
tline = fgetl(fid);
cntr_line = 0;
do_read_para = false;
while ischar(tline)
    cntr_line = cntr_line + 1;
    if contains(tline, 'frameIntervalI')
        do_read_para = true;
    end
    if contains(tline, 'loopIteration')
        NHeader = cntr_line;
        break;
    end
    if do_read_para
        ind = regexp(tline, '",');
        para_name  = tline(2:ind-1);
        para_value = tline(ind+2:end);
        if strcmp(para_value(1), '"')
            try % 'magPID' '"40,,"'
                eval(['para.', para_name, ' = [', para_value(2:end-1), '];'])
            end
        else
            eval(['para.', para_name, ' = ', para_value, ';']);
        end
    end
    tline = fgetl(fid);
end
fclose(fid);

end

