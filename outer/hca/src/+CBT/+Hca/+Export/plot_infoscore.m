function [] = plot_infoscore(barcodeGen )
% plot_infoscore

        figure,
        subplot(2,2,1)
        plot(cellfun(@(x) x.infoscore.mean, barcodeGen),'*')
        xlabel('bar nr')
        title('Mean')
        subplot(2,2,2)
        plot(cellfun(@(x) x.infoscore.std, barcodeGen),'*')
        xlabel('bar nr')
        title('std. deviation')

        subplot(2,2,3)
        plot(cellfun(@(x) x.infoscore.score, barcodeGen),'*')
        xlabel('bar nr')
        title('info score')

        defaultMatFilename ='saveinfoscoreshere';
        [~, matDirpath] = uiputfile('*.txt', 'Save infoscore data as', defaultMatFilename);
        
        
        timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
        matFilename = strcat([timestamp '_' 'means.txt']);
        matFilename2 = strcat([timestamp '_' 'stdev.txt']);
        matFilename3 = strcat([timestamp '_' 'infoscore.txt']);



        if isequal(matDirpath, 0)
            return;
        end

        matFilepath = fullfile(matDirpath, matFilename);

        fileID = fopen(matFilepath,'w');        
        fprintf(fileID,'%6.3e\n',cellfun(@(x) x.infoscore.mean, barcodeGen));
        fclose(fileID);

        matFilepath = fullfile(matDirpath, matFilename2);

        fileID = fopen(matFilepath,'w');        
        fprintf(fileID,'%6.3e\n',cellfun(@(x) x.infoscore.std, barcodeGen));
        fclose(fileID);

        matFilepath = fullfile(matDirpath, matFilename3);

        fileID = fopen(matFilepath,'w');        
        fprintf(fileID,'%6.3e\n',cellfun(@(x) x.infoscore.score, barcodeGen));
        fclose(fileID);

end

