function [bars, bits, names, meanbpnm, pixelWidth_nm, psfSigmaWidth_nm, isLinearTF] = ...
    load_theory_mat(sets,idx)
    %   Args:
    %       - sets -settings containg theory file locations
    %       - idx - index of theory file
    %
    %   Returns:
    %       bars - theory barcodes
    %       bits - bitmasks
    %       names - theory barcode names
    %       meanbpnm - bp/nm parameter
    %       pixelWidth_nm - pixel width nm/px
    %       psfSigmaWidth_nm - psf in nm
    %       isLinearTF - whether barcodes linear   


     %         addpath(genpath(sets.theoryFileFold{idx}));
     if ~isfield(sets,'theoryFile')
        error('No theory file folder selected');
     end

     if isfield(sets,'theoryFileFold')
        barcodeData = load(fullfile(sets.theoryFileFold{idx},sets.theoryFile{idx}));
     else
         barcodeData = load(sets.theoryFile{idx});
     end
        
    bits = [];
    try % there can be different kind of input data, so use the appropriate one
        bars = barcodeData.hcaSessionStruct.theoryGen.theoryBarcodes;
        names = barcodeData.hcaSessionStruct.theoryGen.theoryNames;
        setsB = barcodeData.hcaSessionStruct.theoryGen.sets;
    catch
        try
            bars = barcodeData.hcaSessionStruct.theoryBarcodes;
            names = barcodeData.hcaSessionStruct.theoryNames;
            setsB = barcodeData.hcaSessionStruct.sets;
        catch
            try
                bars = barcodeData.theoryGen.theoryBarcodes;
                bits = barcodeData.theoryGen.theoryBitmasks;
                names = barcodeData.theoryGen.theoryNames;
                setsB = barcodeData.theoryGen.sets;        
            catch
                bars = barcodeData.theoryBarcodes;
                bits = barcodeData.theoryBitmasks;
                names = barcodeData.theoryNames;
                setsB = barcodeData.sets; 
            end
        end 
    end

    if isfield(setsB,'meanBpExtNm')
        meanbpnm = setsB.meanBpExtNm;
        
        pixelWidth_nm = setsB.pixelWidthNm;
        psfSigmaWidth_nm = setsB.psfSigmaWidthNm;
    else
        meanbpnm = setsB.meanBpExt_nm;
        pixelWidth_nm = setsB.pixelWidth_nm;
        psfSigmaWidth_nm = setsB.psfSigmaWidth_nm;
    end
    try  % try to see if setting of linear theory function is included
       isLinearTF  =  setsB.isLinearTF;
    catch
       isLinearTF  =  0;
    end
    
%         theoryStructN = cell(1,length(bars));
%         theoryFold = sets.theoryFileFold{idx};
%         precision = sets.theory.precision;
%         % these settings from old file // if all files are txts, read this
        % from filename

end

