function [matFilepath] = export_plots_txt(idx,fitPositions,firstBarcode,barFit, nameExp,nameTheory,matDirpath) 
    % Exports two plots as txt

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    matFilename = strcat([num2str(idx) '_' 'plots.txt']);
    


    if isequal(matDirpath, 0)
        return;
    end
	matFilepath = fullfile(matDirpath, matFilename);

    fileID = fopen(matFilepath,'w');      
    
    % in case barcode is cut into two places, first index is gonna be
    % larger than the last index
    cutOut  = barFit(fitPositions);
    
    m1 = nanmean(firstBarcode);
    s1= nanstd(firstBarcode);

    m2 = nanmean(cutOut);
    s2= nanstd(cutOut);

    firstBarcode = ((firstBarcode-m1)/s1) *s2+m2;
    
    
    fprintf(fileID,'Output txt for curves\n');
    fprintf(fileID,'[idx] index of selected barcode\n');
    fprintf(fileID,'[fitPositions] start and end pixels where first curve fits along the reference\n');
	fprintf(fileID,'[firstBarcode] intensity profile of the first curve\n');
    fprintf(fileID,'[barFit] intensity curve of the reference\n');
    fprintf(fileID,'[nameExp] name of the experiment curve\n');
    fprintf(fileID,'[nameTheory] name of the theory curve\n');

    % -----
    fprintf(fileID,'%s ',num2str(idx));
    fprintf(fileID,'\n');
    fprintf(fileID,'%s %s ',num2str(fitPositions(1)), num2str(fitPositions(end)));
    fprintf(fileID,'\n');
    fprintf(fileID,'%5d ',firstBarcode);
    fprintf(fileID,'\n');
   

            
    fprintf(fileID,'%5d ',cutOut);
    fprintf(fileID,'\n');
    fprintf(fileID,strcat([nameExp '\n ']));
    fprintf(fileID,strcat([nameTheory '\n ']));

    fclose(fileID);
%     save(matFilepath, 'hcaSessionStruct');
%     
    fprintf('Saved plot data to ''%s''\n', matFilepath);
    
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    matFilename = strcat([num2str(idx) '_' timestamp 'plots_1.txt']);
    matFilepath = fullfile(matDirpath, matFilename);
    fileID = fopen(matFilepath,'w');   
    fprintf(fileID,'%5d\n ',firstBarcode);
    fclose(fileID);
    
   % timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    matFilename = strcat([num2str(idx) '_' timestamp 'plots_2.txt']);
    matFilepath = fullfile(matDirpath, matFilename);
    fileID = fopen(matFilepath,'w');   
    fprintf(fileID,'%5d\n ',cutOut);
    fclose(fileID);
    
    
end
