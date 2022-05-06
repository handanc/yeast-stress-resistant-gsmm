function toptable = readExpressions(experiment)

    %% Read expression toptable
    % Specify range and delimiter
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines = [2, Inf];
    opts.Delimiter = " ";
    % Specify column names and types
    opts.VariableNames = ["Gene", "logFC", "AveExpr", "t", "PValue", "adjPVal", "B"];
    opts.VariableTypes = ["char", "double", "double", "double", "double", "double", "double"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    % Specify variable properties
    opts = setvaropts(opts, "Gene", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "Gene", "EmptyFieldRule", "auto");
    % Read table
    name = sprintf('differential_toptables\\%s_toptable.tsv', experiment);

    toptable = readtable(name, opts);

end