function [myStruct, myNames] = read_default_sets(setstxt)
    % Puts data from setstxt directly to structure

    if nargin < 2
        default = 1;
    end

    if nargin < 3
        depth = 4;
    end

    myStruct = struct();

    if default
        mFilePath = mfilename('fullpath');

        for i=1:depth
            [mFilePath, ~] = fileparts(mFilePath);
        end
        setstxtFull  = fullfile(mFilePath,'files',setstxt);
    end

    try
        setsTable  = readtable(setstxtFull,'Format','%s%s%s','VariableNamingRule','preserve');
    catch
        setsTable  = readtable(setstxt,'Format','%s%s%s','VariableNamingRule','preserve');
    end
    
    cellNames = {};
    myNames = cell(1,size(setsTable,1));
    
    for i = 1:size(setsTable,1)
        validFieldName = split(matlab.lang.makeValidName(setsTable.Var3{i}),'_');
        number = str2num(setsTable.Var1{i});
        myNames{i} = setsTable.Var2{i};
        if ~isnan(number)
            val = number;
        else
            val = strtrim(setsTable.Var1{i});
        end
        if length(validFieldName)>=3 && ~isempty(str2num(validFieldName{2}))
            cellNames = [cellNames, validFieldName{1}];
            myStruct.(validFieldName{1}){str2num(validFieldName{2})} = val; % Example value assignment
        else
            if length(validFieldName)==3
               myStruct.(validFieldName{1}).(validFieldName{2}).(validFieldName{3}) = val; % Example value assignment
            else
                if length(validFieldName)==2
                    myStruct.(validFieldName{1}).(validFieldName{2}) = val; % Example value assignment
                else
                    myStruct.(validFieldName{1}) = val; % Example value assignment
                end
            end
        end
    end
    uniqueNames = unique(cellNames);
    for i=1:length(uniqueNames)
        myStruct.(uniqueNames{i}) = cell2mat(myStruct.(uniqueNames{i}));
    end

    
end