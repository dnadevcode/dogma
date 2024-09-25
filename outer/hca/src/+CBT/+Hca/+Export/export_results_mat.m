function [matFilepath] = export_results_mat(hcaSessionStruct, resultKey)
    % export results
    
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    defaultMatFilename = strcat([resultKey '_' sprintf('%s_%s', timestamp) 'session.mat']);
    [matFilename, matDirpath] = uiputfile('*.mat', 'Save HCA session data as', defaultMatFilename);
    if isequal(matDirpath, 0)
        return;
    end
    matFilepath = fullfile(matDirpath, matFilename);
    save(matFilepath, 'hcaSessionStruct','-v7.3');
    
    fprintf('Saved HCA data ''%s'' to ''%s''\n', resultKey, matFilepath);
end
