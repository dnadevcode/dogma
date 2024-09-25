function [matFilepath] = export_sets_txt(sets, resultKey) 
    % This export the txt file with all the sets files

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    defaultMatFilename = strcat([resultKey '_' sprintf('%s_%s', timestamp) 'session.txt']);
    [matFilename, matDirpath] = uiputfile('*.txt', 'Save chosen settings data as', defaultMatFilename);
    if isequal(matDirpath, 0)
        return;
    end
	matFilepath = fullfile(matDirpath, matFilename);

    fileID = fopen(matFilepath,'w');
    fprintf(fileID,'These are the settings chosen by the user in a HCA session\n\n');    
        
%    fprintf(fileID,'%5d %5d %5d %5d\n',A)
     
    fprintf(fileID,'timeFramesNr %5d\n',sets.timeFramesNr);
    fprintf(fileID,'alignMethod %5d\n',sets.alignMethod);
    fprintf(fileID,'general prestretchMethod %5d\n',sets.prestretchMethod);

    fprintf(fileID,'\n Filter settings\n\n');

    fprintf(fileID,'filter %5d\n',sets.filterSettings.filter);
    fprintf(fileID,'timeFramesNr %5d\n',sets.filterSettings.timeFramesNr);
    fprintf(fileID,'filterMethod %5d\n',sets.filterSettings.filterMethod);
    fprintf(fileID,'filterSize %5d\n',sets.filterSettings.filterSize);

	fprintf(fileID,'\n Barcode consensus settings\n\n');

    
    fprintf(fileID,'barcodeClusterLimit %5d\n',sets.barcodeConsensusSettings.barcodeClusterLimit);
    fprintf(fileID,'barcodeNormalization %s\n',sets.barcodeConsensusSettings.barcodeNormalization);
    fprintf(fileID,'prestretchPixelWidth_nm %5d\n',sets.barcodeConsensusSettings.prestretchPixelWidth_nm);
    fprintf(fileID,'psfSigmaWidth_nm %5d\n',sets.barcodeConsensusSettings.psfSigmaWidth_nm);
    fprintf(fileID,'deltaCut %5d\n',sets.barcodeConsensusSettings.deltaCut);
    fprintf(fileID,'prestretchUntrustedEdgeLenUnrounded_pixels %5d\n',sets.barcodeConsensusSettings.prestretchUntrustedEdgeLenUnrounded_pixels);
    fprintf(fileID,strcat(['stretchFactors ',num2str(sets.barcodeConsensusSettings.stretchFactors) '\n']));
    fprintf(fileID,strcat(['barcodes in consensus ',num2str(sets.barcodeConsensusSettings.barcodesInConsensus) '\n']));
    fprintf(fileID,strcat(['filtered barcodes in consensus ',num2str(sets.filterSettings.barcodesInConsensus) '\n']));

    fprintf(fileID,'\n Barcode generation settings\n\n');
    fprintf(fileID,'meanBpExt_nm %5d\n',sets.barcodeGenSettings.meanBpExt_nm);
    fprintf(fileID,'pixelWidth_nm %5d\n',sets.barcodeGenSettings.pixelWidth_nm);
    fprintf(fileID,'concNetropsin_molar %5d\n',sets.barcodeGenSettings.concNetropsin_molar);
    fprintf(fileID,'concYOYO1_molar %5d\n',sets.barcodeGenSettings.concYOYO1_molar);
        fprintf(fileID,'concDNA %5d\n',sets.barcodeGenSettings.concDNA);

   % fprintf(fileId,'isLinearTF %5d\n',sets.barcodeGenSettings.isLinearTF)
    fprintf(fileID,'deltaCut %5d\n',sets.barcodeGenSettings.deltaCut);
    fprintf(fileID,'widthSigmasFromMean %5d\n',sets.barcodeGenSettings.widthSigmasFromMean);
  %  fprintf(fileID,'yoyo1BindingConstant %5d\n',sets.barcodeGenSettings.yoyo1BindingConstant)
    fprintf(fileID,'computeFreeConcentrations %5d\n',sets.barcodeGenSettings.computeFreeConcentrations);

    fprintf(fileID,'\n Statistics\n\n');
    fprintf(fileID,'pvaluethresh %5d\n',sets.pvaluethresh);

%     save(matFilepath, 'hcaSessionStruct');
%     
    fprintf('Saved chosen settings data ''%s'' to ''%s''\n', resultKey, matFilepath);
end
