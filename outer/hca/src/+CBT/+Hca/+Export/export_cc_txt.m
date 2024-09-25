function [matFilepath] = export_cc_txt(ccvals, lengths,  theoryNames, name,idx, matDirpath) 
    % Exports two plots as txt

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    matFilename = strcat([num2str(idx) '_' timestamp '_' 'ccvals.txt']);
    


    if isequal(matDirpath, 0)
        return;
    end
	matFilepath = fullfile(matDirpath, matFilename);

    fileID = fopen(matFilepath,'w');        

    fprintf(fileID,strcat([name '\n'] ));

    % -----
    fprintf(fileID,'%4.3f, ',ccvals);
    
    fprintf(fileID,'\n');
    for i=1:length(theoryNames)
        fprintf(fileID,strcat([theoryNames{i} '\n'] ));
    end
    
    fprintf(fileID,'%5d, ',lengths);


%     save(matFilepath, 'hcaSessionStruct');
%     
    fprintf('Saved ccvals data to ''%s''\n', matFilepath);
end
