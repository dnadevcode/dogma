function [matFilepath] = export_pvals_txt(pValueResults,pvals,matDirpath) 
    % Exports two plots as txt

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    matFilename = strcat([timestamp '_' 'pvals.txt']);
    %matFilename2 = strcat([timestamp '_' 'rsq.txt']);



    if isequal(matDirpath, 0)
        return;
    end
	matFilepath = fullfile(matDirpath, matFilename);

    fileID = fopen(matFilepath,'w');        


    % -----
    fprintf(fileID,'P-values\n');
    fprintf(fileID,'%6.3e\n',pvals);
    fclose(fileID);
%     save(matFilepath, 'hcaSessionStruct');
%     
    fprintf('Saved pvals data to ''%s''\n', matFilepath);

%     matFilepath = fullfile(matDirpath2, matFilename2);
% 
%     fileID = fopen(matFilepath,'w'); 
%     fprintf(fileID,'R-squared\n');
%     fprintf(fileID,'%6.3e\n',pValueResults.rsq);
%     fclose(fileID);



end
