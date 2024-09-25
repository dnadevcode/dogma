
% load data
idxRun = 1;
oS = resRun{idxRun}.oS;
barcodeGen = bG{idxRun};
minOverlap = 300;

% calculate p-values. TODO: incorporate local similarities (mpmaxlenbased)
% for more accurate results
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS, minOverlap, 0.08,0.01);


% [sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBasedC{idxRun}{1},minLen,thresCC);

curSink = [];
curSource = [];
for i=1:length(sortedVals)
    [curSink(i),curSource(i)] = ind2sub(size(oS),sortedIds(i));
end
sinksource = [curSink; curSource]';


ii = 1;
import Core.comparisonsres_to_struct; 
[rezStruct2] = comparisonsres_to_struct(barcodeGen, comparisonStructC{idxRun},thrS{idxRun}.length);
lenThry = thrS{idxRun}.length;
import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
[f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink(ii),curSource(ii),rezStruct2,lenThry);

[f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,36, 73,rezStruct2,lenThry);
[f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,69, 65,rezStruct2,lenThry);


pos1 = comparisonStructC{idxRun}{curSink(ii)}.pos(1);
pos2 = comparisonStructC{idxRun}{curSource(ii)}.pos(1);
if comparisonStructC{idxRun}{curSource(ii)}.or(1)==1
    posDifTheory = pos1-pos2;
else
    posDifTheory = -(pos1-pos2+comparisonStructC{idxRun}{curSink(ii)}.lengthMatch-comparisonStructC{idxRun}{curSource(ii)}.lengthMatch);
end
posDifExp = oS(curSink(ii),curSource(ii)).pB - oS(curSink(ii),curSource(ii)).pA;

% lenThry = 10000;
posShift = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);

posShift




% ix = sum(~isnan(sortedVals));
ix = 303;
sortedIdsT = sortedIds(1:ix);
sortedValsT = sortedVals(1:ix);

%% new bar-set generation/
import Core.barcode_island_output;
[outConsensus, coverage, consensus, islandElts, islandPx,cGen,badData, dataInfo, barcodeIslandsData, barcodeIslands, barStruct,barIslands]=...
    barcode_island_output(sortedVals,sortedIds, oS, barcodeGen(:)',timestamp,0.02,pxDifSetting, [])
    % 
% 
% import Core.create_barset_from_overlapstruct;
% % [barsetGen] = create_barset_from_overlapstruct(sortedVals, sortedIds, oS, 50, 0.02, barcodeGen(:)');
% 
% import Core.create_graph_from_barsetgen;
% create_graph_from_barsetgen(barsetGen,barcodeGen)

%
% plot the graph
import Core.create_overlap_graph; % overlap graph, without any corrections
[finalgraph,Ggraphs] = create_overlap_graph(size(oS),sortedIdsT,length(sortedIdsT),barcodeGen,zeros(1,length(sortedIdsT)),100,oS);
% print('FIGS/FigS4.eps','-depsc','-r300');
import Core.bar_island_out;
[outConsensus, coverage, consensus,islandElts, islandPx,cGen,badData,dataInfo,barcodeIslandsData,barcodeIslands,barStruct,barIslands] = bar_island_out(sortedValsT,sortedIdsT, oS,barcodeGen,timestamp,[],50);




    nonEmpty= find(cellfun(@(x) ~isempty(x),barIslands)); % keep only non-empty
barIslands = barIslands(nonEmpty)
cellfun(@(x) x.curPair,barIslands,'UniformOutput',false)


cellfun(@(x) max(x(:,2))-min(x(:,1)),barcodeIslandsData)

import CBT.Hca.UI.Helper.plot_any_bar;
plot_any_bar(10,bG{idxRun},rezMax,thrS{idxRun},1);
 

%
% [mpMax,finalgraph,Ggraphs,sortedVals,sortedIds,pscores] = bargrouping_funs(oS,bG,0.8)
% toc
ixC = [sub2ind(size(oS),78, 86);...
    sub2ind(size(oS),9,86);...
    sub2ind(size(oS),9,78);...
    sub2ind(size(oS),12,78);...
    ];

vals = arrayfun(@(x) find(sortedIds==x),ixC);

sortedIdsA = sortedIds(vals);
sortedValsA = sortedVals(vals);
import Core.bar_island_out;
[outConsensus, coverage, consensus,islandElts, islandPx,cGen,badData,dataInfo,barcodeIslandsData,barcodeIslands,barStruct] = bar_island_out(sortedValsA,sortedIdsA, oS,barcodeGen,timestamp);
% print('FIGS/Fig2.eps','-depsc','-r300');
