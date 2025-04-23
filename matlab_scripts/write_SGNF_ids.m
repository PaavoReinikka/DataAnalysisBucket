function y = write_SGNF_ids(fulldata, h, TAG)
    data = fulldata(h,:);
    ids = split(data.VarName1, '_');
    filenames = ids(:,1);
    unique_names = unique(filenames);
    
    for i=1:length(unique_names)
        % names of a file
        file = unique_names(i);
        % rowinfo of SGNF features
        info = split(data.VarName1(filenames==file,1),'_');
        % rownumber (only uniques)
        id = unique(str2double(info(:,3)));
        % filename for saving
        fname = sprintf("%s%s_ID.csv", TAG, file);
        % make a table and save to file
        T = table(id);
        writetable(T, fname, 'Delimiter',';','WriteVariableNames',0);

        fprintf("File: %s  -  %d\n", file, length(id))
    end

end