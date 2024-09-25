function [outputArg1,outputArg2] = merge_validation(inputArg1,inputArg2)


[f] = pair_evaluation_with_ground_truth_plot(barcodeGen', oS,2,20,synthStr2{idxRun},lenThry);

[f] = pair_evaluation_with_ground_truth_plot(barcodeGen', oS,20,2,synthStr2{idxRun},lenThry);


bId = barcodeIslands{2}(1:3);
bars = barcodeGen(bId);

oSsub = synthStr2{idxRun}([bId],[bId]);
[~] = bar_island_out(sortedValsGoodSub,sortedIdsGoodSub, oSsub,bars',timestamp);

oSsub = oS([bId],[bId]);
[~] = bar_island_out(sortedValsGoodSub,sortedIdsGoodSub, oSsub,bars',timestamp);


[sortedValsSub, sortedIdsAllSub, pscoresSub,foSSub,overlaplenSub] = sorted_scores(oSsub);

    [passLocalThreshSub] = filter_scores_local(sortedValsSub,mpMaxLenBased,minOverlap,minLen,3);
    [passGlobalThreshSub] = filter_scores_global(foSSub,mpMaxLenBased,overlaplenSub,minLen,3);
    
goodPosSub = passLocalThreshSub.*passGlobalThreshSub;
sortedIdsGoodSub = sortedIdsAllSub(find(goodPosSub));
sortedValsGoodSub =sortedValsSub(find(goodPosSub));

import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(size(oSsub),sortedIdsGoodSub,length(sortedIdsGoodSub),bars,zeros(1,length(sortedIdsGoodSub)),100,oSsub);



[~] = bar_island_out(sortedValsGoodSub,sortedIdsGoodSub, oSsub,bars',timestamp);



figure,plot(bars{1}.rawBarcode)
hold on
plot(bars{2}.rawBarcode+5)
plot(bars{3}.rawBarcode+10)

end

