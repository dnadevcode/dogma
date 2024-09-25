function [] = a_validation()
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

idxRun = 2;
oS = resRun{idxRun}.oS;
barcodeGen = bG{idxRun};
lenThry = theoryStruct{1}{1}.length;
minOverlap = 300;
% 
Nvals = 2000;

[sortedVals, sortedIdsAll, pscores,foS,overlaplen] = sorted_scores(oS);
posShift = calc_pos_dif(oS,sortedIdsAll,synthStr{idxRun},lenThry, Nvals);

[passLocalThresh] = filter_scores_local(sortedVals,mpMaxLenBased,minOverlap,minLen,Nvals);
[passGlobalThresh] = filter_scores_global(foS,mpMaxLenBased,overlaplen,minLen,Nvals);

goodPos = passLocalThresh.*passGlobalThresh;
sortedIdsGood = sortedIdsAll(find(goodPos));
sortedValsGood =sortedVals(find(goodPos));

import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(size(oS),sortedIdsGood,length(sortedIdsGood),barcodeGen,zeros(1,length(sortedIdsGood)),100,oS);


% refinement to make graph circular? do we need to remove weak links (that
% might represent false mappings?)

% plot output
[outConsensus, coverage, consensus, islandElts, islandPx,cGen,badData,dataInfo,barcodeIslandsData,barcodeIslands] = bar_island_out(sortedValsGood,sortedIdsGood, oS,barcodeGen',timestamp);
figure,plot(cellfun(@(x) x.positionalDifference,dataInfo),'x')

% true possitions (as calculated via synthStr2)
[outConsensus2, coverage2, consensus2, islandElts2, islandPx2,cGen,badData,dataInfo,barcodeIslandsData,barcodeIslands] = bar_island_out(sortedValsGood,sortedIdsGood, synthStr2{idxRun},barcodeGen',timestamp);

% how to evaluate correctness for synthetic data?
% check the scaling factors. Should be mean=1 if data is consistent, if
% it's larger or smaller, we might want to re-scale to have = 1 (should
% give more accurate overall length)
figure,histogram(cellfun(@(x) x.bestBarStretch,cGen{1}.comparisonStruct))

%% note: badDataInfo with very large poistionalDifference parameter (N~10000) indicates the left and right edge of the consensus barcode, and it should be 
% circular through this location! One could get an alternative assembly by
% starting from this

%%
w = 400;% minimum overlap length
sF = 0.9:0.025:1.1;

bT{1}.rawBarcode = theoryStruct{idxRun}{1}.rawBarcode;
bT{1}.rawBitmask = ones(1,length(bT{1}.rawBarcode));
bT{1}.rawBitmask = logical(bT{1}.rawBitmask);

wBoostrapping = 2000;
numWindows = floor(length( cGen{1}.rawBarcode)/wBoostrapping);
R = reshape(1:floor(length(cGen{1}.rawBarcode)/numWindows)*numWindows,floor(length(cGen{1}.rawBarcode)/numWindows),[]);
     

f = figure;
tiledlayout(size(R,2),2,'TileSpacing','tight');

for k=1:size(R,2);
    cos{1}.rawBarcode = cGen{1}.rawBarcode(R(:,k));
    cos{1}.rawBitmask = logical(cGen{1}.rawBitmask(R(:,k)));
    
    %  add theory as last barcode
    consensusGen =[cos bT];
    
    % calculate MP scores 
    [oScons{1}] = calc_overlap_mp(consensusGen,sF, w,timestamp);
     
    % [oScons] = calc_scores({cGen},w,sF);
    
    
           
    %% Plot specific pair
    
    import Core.plot_match_simple;
    expidx = 1; % idx of experiment
    thrIdx = 2; % idx of theory
    lenThry = theoryStruct{idxRun}{1}.length;
    import Plot.pair_evaluation_with_ground_truth_plot;
    [f] = pair_evaluation_with_ground_truth_plot(consensusGen, oScons{1},expidx,thrIdx,[],lenThry,f);
end

% [f] = plot_match_simple(consensusGen, oScons{1},expidx,thrIdx);
% [f] = plot_match_simple(cGen, oScons{1},numExpBar+thrIdx,expidx);

%
%% plot comparison for the best pair (this is both for theory vs exp and exp vs exp
% idFold = 1;
% idxPair = 3; % 1 is the pair with best score
%  
% 
% 
% import Core.plot_match_pcc;
% [sortedVals, sortedIds,pscores] = sorted_scores(oScons{idFold});
% 
% 
% [curSink,curSource] = ind2sub(size(oScons{idFold}),sortedIds(idxPair));
% 
% 
% import Core.plot_match_simple;
% % [f] = plot_match_simple(barStruct, oS,curSink,curSource);
% [f] = plot_match_simple(cGen, oScons{idFold},curSink,curSource);



end

