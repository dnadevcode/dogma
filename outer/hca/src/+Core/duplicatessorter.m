function [duplicateInfo,oS] = duplicatessorter(hcaSets)
    %   
    %   Args:
    %       hcaSets - settings
    %
    %   Returns:
    %       duplicateInfo - info about duplicates
    %       oS - resulting oS information


    import CBT.Hca.Import.add_kymographs_fun;
    import CBT.Hca.Core.align_kymos;

    if nargin < 2
        hcaSets.kymosets.filenames = hcaSets.kymofolder;
        hcaSets.kymosets.kymofilefold = cell(1,length(   hcaSets.kymosets.filenames));
       % add kymographs
        [kymoStructs] = add_kymographs_fun(hcaSets);
    end

    sF = 1; % length re-scaling factors. Define via sets

    % keep single timeframe
    import CBT.Hca.Core.edit_kymographs_fun;
    kymoStructs = edit_kymographs_fun(kymoStructs,1);

    [ kymoStructs ] = align_kymos( hcaSets, kymoStructs );

    barcodeGen = cell(1,length(kymoStructs));
    for i=1:length(kymoStructs)
        barcodeGen{i}.rawBarcode = kymoStructs{i}.alignedKymo;
        barcodeGen{i}.rawBitmask = logical(kymoStructs{i}.alignedMask);
    end

    disp(['Running all-to-all pairwise comparison for ', num2str(length(kymoStructs) ), ' barcodes'])

    % pairwise comparison using full overlap PCC
    tic
    import  CBT.Hca.Core.Comparison.compare_pairwise_distance;
    oS = compare_pairwise_distance(barcodeGen,sF, hcaSets.default.minLen);
    timeUsed = toc;

    localScore = [oS(:).sc]; % local overlap score

    %


    %

    lenA = [oS(:).lenA]; % lenA
    lenB = [oS(:).lenB]; % lenB
    overlaplen = [oS(:).overlaplen];


    idx = reshape(1:size(oS,1)*size(oS,2), size(oS,1),size(oS,2)); % from bg_test_1
    idx = tril(idx);

    localScore(idx(idx~=0)) = nan;

    locsDuplicates = find(localScore > hcaSets.default.Cthresh);

    % should we use it as well, or it skips some?
    allOverlaps = overlaplen./max(lenB,lenA);
    
    % allOverlaps(locsDuplicates)

    bar1 = zeros(1,length(locsDuplicates));
    bar2 = zeros(1,length(locsDuplicates));
    names = cell(1,length(locsDuplicates));

   for i=1:length(locsDuplicates)
        [bar1(i), bar2(i)] = ind2sub(size(oS),locsDuplicates(i));
        names{i} = {kymoStructs{bar1(i)}.name; kymoStructs{bar2(i)}.name};
   end

    duplicateInfo.locsDuplicates = locsDuplicates;
    duplicateInfo.bar1 = bar1;
    duplicateInfo.bar2 = bar2;
    duplicateInfo.names = names;

    disp(['Found [', num2str(duplicateInfo.locsDuplicates ), '] duplicates'])

    save('duplicate.mat',duplicateInfo)


%     partialScore(idx(idx~=0)) = nan;

    % plot duplicates:
    %     ix = 1;
    %     import Core.plot_match_pcc;
    %     [f] = plot_match_pcc(barcodeGen, oS,bar1(ix),bar2(ix));



end

