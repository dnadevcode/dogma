function [theoryStruct] = theory_mat_to_struct(theoryGen)
    % convert theory mat output to structure with proper fields

    bars = theoryGen.theoryBarcodes;
    bits = theoryGen.theoryBitmasks;

    names = theoryGen.theoryNames;
    setsB = theoryGen.sets;
 
    meanbpnm = setsB.meanBpExt_nm;
    pixelWidth_nm = setsB.pixelWidth_nm;
    psfSigmaWidth_nm = setsB.psfSigmaWidth_nm;
    try  % try to see if setting of linear theory function is included
       isLinearTF  =  setsB.isLinearTF;
    catch
       isLinearTF  =  0;
    end
    
    
    theoryStruct = cell2struct([bars;...
    bits;arrayfun(@(x) isLinearTF,1:length(bars),'un',false);...
    names;...
    cellfun(@(x) length(x),bars,'un',false);...
    arrayfun(@(x) meanbpnm,1:length(bars),'un',false);...
    arrayfun(@(x) pixelWidth_nm,1:length(bars),'un',false);...
    arrayfun(@(x) psfSigmaWidth_nm,1:length(bars),'un',false)]',...
    {'rawBarcode','rawBitmask', 'isLinearTF','name','length','meanBpExt_nm','pixelWidth_nm','psfSigmaWidth_nm'},2);
     

end

