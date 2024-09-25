function [barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGen,barcodeIslandsData, barStruct,barIslands]=...
    barcode_island_output(sortedScores,sortedIds, oS, bG,timestamp,scDiffSetting,pxDifSetting, tsVisual,plotFig, N)
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

    if nargin < 8
        tsVisual = [];
    end

    if nargin < 9 || isempty(plotFig)
        plotFig = 1;
    end

    if nargin < 10
        N = 2; % remove only the ones containing two elements
    end

    barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),bG,'un',false);...
        cellfun(@(x) x.rawBitmask,bG,'un',false)]',{'rawBarcode','rawBitmask'},2);
    

    % Most important algorithm, creating barsets
    
    %% new bar-set generation/
    import Core.create_barset_from_overlapstruct;
    [barsetGen] = create_barset_from_overlapstruct(sortedScores, sortedIds, oS, pxDifSetting, scDiffSetting , barStruct);

            


        
    % find groups with more than N members
    nonEmpty= find(cellfun(@(x) length(x.members)>N,    barsetGen.barcodeIslands)); % keep only non-empty
    %     ix =4;
%     import Plot.plot_single_island;
%     plot_single_island(barsetGen.barcodeIslandsData{nonEmpty(ix)},barsetGen.barcodeIslands{nonEmpty(ix)}.members,barStruct)

% 
%     barcodeIslandsData = barcodeIslandsData(nonEmpty);
%     barcodeIslands = barcodeIslands(nonEmpty);
%     barislandLength = cellfun(@(x) max(x(:,2))-min(x(:,1)),barcodeIslandsData)

    barIslands = arrayfun(@(x)  barsetGen.barcodeIslands{x}.members,nonEmpty,'un',false);
    barcodeIslandsData =  arrayfun(@(x)  barsetGen.barcodeIslandsData{x},nonEmpty,'un',false);
%     % plot group quick

   % vs. consensus positions
    import Plot.islandsPosStruct;
    posStruct = islandsPosStruct(barcodeIslandsData,barIslands);


%     import Plot.plot_best_pos;
%     fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData))),1);%1/bpPx*10^6

    
    [outConsensus, coverage, consensus,islandElts, islandPx] = ...
        gen_assembly(bG(cellfun(@(x) x.barid,posStruct))',posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData))),timestamp,0);
    
    if plotFig
        import Plot.plot_island;
        plot_island(posStruct,islandElts, bG(cellfun(@(x) x.barid,posStruct))',1,tsVisual);
    end

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