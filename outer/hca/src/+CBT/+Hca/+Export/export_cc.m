function [matFilepath,timestamp] = export_cc(T, matDirpath, timestamp, writemode) 
    % Exports table as txt

    if nargin < 3
        timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
        writemode = 'overwrite';
    end
    matFilename = strcat([ 'table_' timestamp  '.txt']);
    matFilename2 = strcat([ 'table_' timestamp  '.dat']);


    if isequal(matDirpath, 0)
        return;
    end
	matFilepath = strcat([matDirpath, matFilename]);
    writetable(T,matFilepath,'Delimiter','\t','WriteMode',writemode)  
    matFilepath = strcat([matDirpath, matFilename2]);

    writetable(T,matFilepath,'WriteRowNames',true,'WriteMode',writemode)  
    fprintf('Saved table to ''%s''\n', matFilepath);
end
