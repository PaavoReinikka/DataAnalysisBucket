function allModes = importfile(filename)
    %IMPORTFILE Import data from a text file
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 49);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ";";
    
    % Specify column names and types
    opts.VariableNames = ["VarName1", "VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "aSYN", "aSYN1", "aSYN2", "aSYN3", "aSYN4", "aSYN5", "aSYN6", "aSYN7", "aSYN8", "aSYN9", "comb", "comb1", "comb2", "comb3", "comb4", "comb5", "comb6", "comb7", "comb8", "comb9", "INFg", "INFg1", "INFg2", "INFg3", "INFg4", "INFg5", "INFg6", "INFg7", "INFg8", "INFg9", "UT", "UT1", "UT2", "UT3", "UT4", "UT5", "UT6", "UT7", "UT8", "UT9"];
    opts.VariableTypes = ["string",   "string",   "string",   "double",   "double",   "double",   "double",  "categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["VarName1", "VarName2", "VarName3", "VarName8"], "EmptyFieldRule", "auto");
    %opts = setvaropts(opts, "VarName9", "TrimNonNumeric", true);
    %opts = setvaropts(opts, "VarName9", "ThousandsSeparator", ",");
    
    % Import the data
    allModes = readtable(filename, opts);

end