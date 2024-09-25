function [outConsensus, coverage, consensus, islandElts, islandPx,cGen,badData, dataInfo, barcodeIslandsData, barcodeIslands, barStruct,barIslands]=...
    bar_island_out(sortedScores,sortedIds, oS, bG,timestamp,tsVisual,pxDifSetting)
    %   Args:

    %   sortedScores - sorted scores
    %   sortedIds - sorted ids corresponding to scores
    %   oS - overlap structure
    %   bG - bargroup that we're analysing
    %   timestamp - timestamp for current experiment
    %   tsVisual - visual output
    %   pxDifSetting - how many pixels should maps be in order to merge
    %   them
    %

    %   Returns:
    %
    %       outConsensus - output consensus
    %       coverage -  output coverage

    % TODO: simplify output. Possibly define scoreRegistry,
    % scoreRegistryIndex, outputTree as for the previous method

    if nargin < 6
        tsVisual = [];
    end

    barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),bG,'un',false);...
        cellfun(@(x) x.rawBitmask,bG,'un',false)]',{'rawBarcode','rawBitmask'},2);
    

    % Most important algorithm, creating barsets
    import Core.create_barset_os;
    [barcodeIslands,barcodeIslandsData, badData,dataInfo, barIslands] = ...
        create_barset_os(sortedScores,sortedIds, barStruct,oS,[],pxDifSetting);

%     import Core.create_barset_from_overlapstruct;
%     [barcodeIslands, barcodeIslandsData, badData, dataInfo, barIslands] = create_barset_from_overlapstruct(sortedScores, sortedIds, oS, pxDifSetting)

    
    nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands)); % keep only non-empty
    barcodeIslandsData = barcodeIslandsData(nonEmpty);
    barcodeIslands = barcodeIslands(nonEmpty);
    barislandLength = cellfun(@(x) max(x(:,2))-min(x(:,1)),barcodeIslandsData)


    % plot group quick
%     sourceGroup = 1;
%     plot_temp(barcodeIslandsData{sourceGroup}, barcodeIslands{sourceGroup},barStruct )

    
    % assembly w/o reference / bars are shuffled..
    nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
    import Plot.islandsPosStruct;
    posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
    % import Plot.plot_best_pos;
    % fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6

    
    [outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(bG(cellfun(@(x) x.barid,posStruct))',posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp,0);
    
    
    import Plot.plot_island;
    plot_island(posStruct,islandElts, bG(cellfun(@(x) x.barid,posStruct))',1,tsVisual);

    allIdx = cellfun(@(x) x.barid,posStruct);

    % put consensus barcodes into a structure
    cGen = cell(1,length(consensus));
    for i=1:length(consensus)
        cGen{i}.rawBarcode = consensus{i};
        cGen{i}.rawBitmask = ~isnan(consensus{i});
        cGen{i}.comparisonStruct = posStruct(islandElts{i}); % is this correct, or are these barid in posStruct?
        cGen{i}.idx = allIdx(islandElts{i});
    end

%% todo: create a file with position/orientation of each barcode wrt consensus


end