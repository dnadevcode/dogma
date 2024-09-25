% load mat

%%
dataFold = {'/export/scratch/albertas/data_temp/bargrouping_selected/'};
dataSets = {'ecoli_S2/','ecoli_S5/','ecoli_dEC','synth','ecoli_test/','ecoli_test2/'}; % different available datasets
fastas = {'018_final_polish.fasta','DA32087.fasta'};

% 2-nd probably not 1
theoryIdxs = {1, 1, 2, nan, 2 , 2}; % known theory indexes
dataSetIdx = 3; % 5,... test datasets (part of the data)


nmbp = mean(cellfun(@(x)   x{1}.nmBpidFold,kymoStructs));


%%
sets.minOverlap = 300;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, 0.08,0.01); %localCCStruct


sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.04;
sets.pxDifSetting = 50;

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [])

import Core.create_graph_from_barsetgen;
 [Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen);

meanSF = mean(barcodeIslandsData{1}(:,end));

approxNMBP = meanSF*nmbp

% validation against theory
import Validation.barislands_vs_theory_validation_all

barislands_vs_theory_validation_all



% nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty
% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty


%% Need to compare barcodes at specific NM/BP
% sets.comparisonMethod = 'mass_pcc';
% [comparisonStructC{1},rezMax,bestBarStretch] = compare_to_t(barcodeGen,theoryStruct,0.9:0.01:1.1,sets);


%%
idxPair = 50;
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));

% lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
import Plot.pair_evaluation_plot;
[f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,comparisonStructInd([curSink curSource]),sum(thrInd.length));
%% plot
% print('FIGS/Fig1.eps','-depsc','-r300');

%%
N = 1000;
passThreshPair = zeros(1,N);
posShift = zeros(1,N);
for idxPair=1:N
    
    [curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
%     if ismember(curSink,passthreshC{idxRun})&&ismember(curSource,passthreshC{idxRun})
%         passThreshPair(idxPair) = 1;
%     end
    pos1 = comparisonStructInd{curSink}.pos(1);
    pos2 = comparisonStructInd{curSource}.pos(1);
    if comparisonStructInd{curSource}.or(1)==1
        posDifTheory = pos1-pos2;
    else
        posDifTheory = -(pos1-pos2+comparisonStructInd{curSink}.lengthMatch-comparisonStructInd{curSource}.lengthMatch);
    end
    posDifExp = oS(curSink,curSource).pB - oS(curSink,curSource).pA;


    
    posShift(idxPair) = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);
end

%   passthreshC{idxRun}  these passed thresh
figure,tiledlayout(2,1)
nexttile
plot(posShift);hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
% lgn = legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})
nexttile
plot(posShift(1:length(barsPassThreshMP)+50));hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
plot(find(~passThreshPair(1:length(barsPassThreshMP)+50)),posShift(find(~passThreshPair(1:length(barsPassThreshMP)+50))),'redx')
lgn = legend({'Distance between estimated GT and pair overlap','Threshold for good scores','Not accurate map against theory'});
lgn.Location ='southoutside';
% print('FIGS/FigS3.eps','-depsc','-r300');



      
import Validation.calc_pos_dif;
posShift = calc_pos_dif(oS,sortedIdsGood(1:N),comparisonStructInd,lenThry, N);

import Plot.position_evaluation_plot_stouffer;
[fp] = position_evaluation_plot_stouffer(sortedValsGood,localScore,posShift)

%         lenThry = 