function [  ] = pregenerate_cc_db( path, newNmBp,barnr, lenBarcodes)

    if nargin < 2
        path = 'allbarcodes_2018-09-05_11_31_46_session.mat';
        newNmBp = 0.23;
    end

    if nargin < 4 
        barnr = 24;
        lenBarcodes = 300;
    end
    % load theory
    import CBT.Hca.Import.load_thry;
    [theoryStruct, sets ] = load_thry(path,newNmBp);

    sets.barcodeConsensusSettings.aborted = 1;
    % no filtering
    sets.filterSettings.filter = 0;

    % stretch factors, could be passed as an input
    sets.barcodeConsensusSettings.stretchFactors  = [1];
    sets.filterSettings.filterSize = 2.3;

    import CBT.Hca.Core.Comparison.compare_to_thry;
    [ maxVals ] = compare_to_thry( theoryStruct, barnr, lenBarcodes, sets );

    rS = 'ccMaxChrom.txt';
    fid = fopen(rS,'w');
    for i = 1:length(maxVals)
        fprintf(fid, '%5d ', maxVals(i));
        fprintf(fid, '\n');
    end
    fclose(fid);


end

