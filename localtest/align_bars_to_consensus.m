
idx = 1;

consensusbar = consensus{idx};

% align all bars to consensusbar


bars = barcodeGen(barIslands{idx}); 

barx = barIslands{idx};

%todo: unique sf for each bar

sets.comparisonMethod = 'mass_pcc';

%% Can loop this to re-align barcodes:

tS.rawBarcode = consensusbar(~isnan(consensusbar ));
tS.rawBitmask = ones(1,length(tS.rawBarcode)); 
tS.length = sum(~tS.rawBitmask);
tS.isLinearTF = 1;
tS.name = 'test';
% tS.rawBarcode(isnan(tS.rawBarcode))

sF = 0.9:0.01:1.1;

% compare to theory: either pcc or mp
[comparisonStruct,rezMax,~] = compare_to_t(bars,tS,sF,sets);

import Core.comparisonsres_to_struct; 
[oSreal] = comparisonsres_to_struct(bars, comparisonStruct,tS.length);

% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oSreal,nan, 0.08,0.01); %localCCStruct

sortedVals = sortedVals(1:length(localScore));
sortedVals(100:end) = nan; % keep less scores to show two islands
sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.02;
sets.pxDifSetting = 20; %20

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oSreal, bars,timestamp,sets.scDiffSetting,sets.pxDifSetting, [])

