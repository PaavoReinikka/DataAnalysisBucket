function y = write_new_data(fulldata, h, TAG)
    data = fulldata(h,:);
    ids = split(data.VarName1, '_');
    filenames = ids(:,1);
    unique_names = unique(filenames);
    
    for i=1:length(unique_names)
        file = unique_names(i);
        fname = sprintf("%s%s.csv", TAG, file);
        T = data(filenames==file,2:end);
        writetable(T, fname, 'Delimiter',';');
    end

end