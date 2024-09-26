function [] = final_figure_4(sets, resRun, bG,synthStr,theoryStruct)
% figure 4

% load('C:\Users\Lenovo\postdoc\DATA\TEMP\5_2023-11-16_20_51_27resRun.mat');
alphaNu1 = 0.09;
alphaNu2 = 0.085;
alphaN1 = 0.42;
alphaN2 = 0.004;

pthresh = 0.05;
nB = 200;
sets.minOverlap = 300;
sets.scDiffSetting = 0.05;
sets.pxDifSetting = 50;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

ixtest = 8; %setsToRun(5) % load the synthetic dataset with 3 islands




% maybe put into sep function
oS = resRun{ixtest}.oS;
barcodeGen = bG{ixtest}';

oS = oS(1:nB,1:nB);
barcodeGen = barcodeGen(1:nB);

% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu1, pthresh,[],alphaN1,alphaNu2,alphaN2); %localCCStruct

sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood, sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],1,2)
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen,[],0);


%%
N = length(sortedValsGood);

import Validation.calc_distance_scaled; % need to properly scale things when calculating distances for validation
posShift = calc_distance_scaled(oS,sortedIdsGood(1:N), synthStr{ixtest},theoryStruct{ixtest}{1}.length, N);

N = min(1000,length(sortedIdsBad));


import Validation.calc_pos_dif;
posShiftBad = calc_distance_scaled(oS,sortedIdsBad(1:N), synthStr{ixtest},theoryStruct{ixtest}{1}.length, N);


% calc_distance_scaled

% import Validation.calc_pos_dif;
% posShift = calc_pos_dif(oS,sortedIdsGood(1:N), synthStr{ixtest},theoryStruct{ixtest}{1}.length, N);


% 
% import Validation.calc_pos_dif;
% posShiftBad = calc_pos_dif(oS,sortedIdsBad(1:N), synthStr{ixtest},theoryStruct{ixtest}{1}.length, N);
%%
f=figure('Position', [10 10 514 500]),g=tiledlayout(3,2,'TileSpacing','tight','Padding','none')



idxPair = 1;
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
import Plot.pair_evaluation_with_ground_truth_simple_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,curSink,curSource,[],theoryStruct{ixtest}{1}.length,g);


import Plot.position_evaluation_stouffer;
[f,g] = position_evaluation_stouffer(sortedValsGood,sortedValsBad, localScore,posShift,posShiftBad,f,g);



% import Plot.pair_evaluation_ground_truth; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
% pair_evaluation_ground_truth(barcodeGen, oS,curSink,curSource,synthStr{ixtest},theoryStruct{ixtest}{1}.length,g);


print('FIGS/Fig4.eps','-depsc','-r500');

[sum(posShift>30)/length(posShift) mean(posShift(posShift<30)) std(posShift(posShift<30))]

f=figure('Position', [10 10 900 200]);g=tiledlayout(1,2,'TileSpacing','tight','Padding','none');

import Plot.pair_evaluation_ground_truth; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_ground_truth(barcodeGen, oS,curSink,curSource,synthStr{ixtest},theoryStruct{ixtest}{1}.length,g);
print('FIGS/Fig4S.eps','-depsc','-r500');

% 
% import Plot.position_evaluation_stouffer;
% [fp] = position_evaluation_stouffer(sortedValsGood,sortedValsBad, localScore,posShiftBad)
% 
% 
% 
% import Plot.position_evaluation_plot_stouffer;
% [fp] = position_evaluation_plot_stouffer(sortedIdsBad,localScore,posShift,posShiftBad)
% 
% 
% 
% 
% f = figure
% g =tiledlayout(4,2,'TileSpacing','compact')
% 
% idxPair = 1;
% [curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
% 
% 
% import Plot.pair_evaluation_with_ground_truth_simple_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
% pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,curSink,curSource,synthStr{ixtest},theoryStruct{ixtest}{1}.length,g);
% 
% 
% idxPair = 100;
% [curSink,curSource] = ind2sub(size(oS),sortedIdsBad(idxPair));
% 
% 
% import Plot.pair_evaluation_with_ground_truth_simple_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
% pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,curSink,curSource,synthStr{ixtest},theoryStruct{ixtest}{1}.length,g);

% 
% import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
% pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,synthStr2{ixtest},theoryStruct{ixtest}{1}.length,g);

% import Plot.pair_evaluation_plot;
% [f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,synthStr{ixtest}([curSink curSource]),theoryStruct{ixtest}{1}.length);
% %

end

