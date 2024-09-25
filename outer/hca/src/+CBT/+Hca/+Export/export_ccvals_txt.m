function [matFilepath] = export_ccvals_txt(ccvals,matDirpath) 
    % Exports two plots as txt

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    matFilename = strcat([timestamp '_' 'ccvals.txt']);
    


    if isequal(matDirpath, 0)
        return;
    end
	matFilepath = fullfile(matDirpath, matFilename);

    fileID = fopen(matFilepath,'w');        


    % -----
    fprintf(fileID,'%4.3f, ',ccvals);

%     save(matFilepath, 'hcaSessionStruct');
%     
    fprintf('Saved ccvals data to ''%s''\n', matFilepath);
end
