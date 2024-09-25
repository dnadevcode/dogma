function [ theorBar,theorBit, expBar, expBit] = load_theory_and_stretch_ex(ii, theoryStruct, comparisonStruct,barcodeGen )
    % load_theory_and_stretch_ex
    
    % load theory file
    fileID = fopen(theoryStruct{comparisonStruct{ii}.idx}.filename,'r');
    formatSpec = '%f';
    theorBar = fscanf(fileID,formatSpec)';
    fclose(fileID);
    
    % theory length
    thrLen = theoryStruct{comparisonStruct{ii}.idx}.length;
    
    % bitmask. In case of linear barcode, would like to modify this
    theorBit = ones(1,thrLen);
    
    % load either theory barcode or the consensus barcode
    try
        expBar = barcodeGen{ii}.rawBarcode;
        expBit = barcodeGen{ii}.rawBitmask;
    catch
        expBar = consensusStruct.rawBarcode;
        expBit = consensusStruct.rawBitmask;  
    end
    expLen = length(expBar);

    % interpolate to the length which gave best CC value
    expBar = interp1(expBar, linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch));
    expBit = expBit(round(linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch)));
    expBar(~expBit)= nan;
    


end

