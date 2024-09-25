function [] = ol_validation(synthStr,theoryStruct, bG, sF, minOverlap)


timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% sF = 0.95:0.01:1.05; % length re-scaling factor
sF = 0.9:0.01:1.1; % length re-scaling factor

minOverlap = 300; % minimum overlap
% minOverlapContig = 600;
% minLen = 300;

plotfigs=0;

N = 1000;

%% mpMaxLenBased - possibly make a number of runs and get
% nmPerPx = 110;
% nmbp = 0.2 % length re-scale factor
% minL = 300;
% tL = theoryStruct{1}{1}.length;
% minLen =[minL:50:3000]; % min to max
% [MP,mpMaxLenBased,theoryStructRev,MPI,~] = bargrouping_minimum_length([],nmPerPx,nmbp,1,30,minLen, tL*nmPerPx/nmbp);

import Core.plot_match_pcc;
import Plot.pair_evaluation_plot;
import Plot.pair_evaluation_with_ground_truth_plot;
import Validation.calc_pos_dif;
import Nullmodel.filter_scores_local; % all scores up to a value, since  minOverlap is fixed   
import Nullmodel.filter_scores_global; % all scores up to a value, since  minOverlap is fixed   

resRun = cell(1,length(bG));
for idxRun = 1:length(bG)
    barcodeGen = bG{idxRun};
    lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(barcodeGen',sF, minOverlap, timestamp);
    toc

% just  PCC
%     import  CBT.Hca.Core.Comparison.compare_pairwise_distance;
%     oS = compare_pairwise_distance(barcodeGen,sF, hcaSets.default.minLen);

    resRun{idxRun}.oS = oS;
end
save(['synth','_',timestamp,'resRun.mat'] ,'resRun','bG','synthStr2','synthStr','theoryStruct','minOverlap','sF','setsGen')

  



 c = parcluster;
% c.AdditionalProperties.AccountName = 'snic2021-5-132';
c.AdditionalProperties.WallTime = '6:00:00';
% j =batch(c,@local_pipeline_mp,1,{ix, iy, dirName, depth, 0, sF, thryFiles},'Pool',30);



%%
idxRun = 3;
oS = resRun{idxRun}.oS;
barcodeGen = bG{idxRun}';
% [sortedVals, sortedIdsAll, pscores,foS,overlaplen,partialScore,partialLength,~,~,scoreCombined] = sorted_scores(oS);
  

% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, 0.08,0.01); %localCCStruct


sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.02;
sets.pxDifSetting = 20;

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [])
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty



    import Core.create_graph_from_barsetgen;
    create_graph_from_barsetgen(barsetGen,barcodeGen)


    
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS,minOverlap, 0.07,0.01);



import Nullmodel.full_dist;

    % plot1: take ground trtuh comparisonStruct

    if plotfigs
%         figure,histogram(foS)

        idxPair = 4;
        [curSink, curSource] = ind2sub(size(oS),sortedIdsGood(idxPair));

%         lenThry = theoryStruct{idxRun}{1}.length;
%         [f] = pair_evaluation_plot(barcodeGen', oS,curSink,curSource,synthStr{idxRun}([curSink curSource]),lenThry);
        
        lenThry = theoryStruct{idxRun}{1}.length;
        [f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,synthStr2{idxRun},lenThry);
    print('FIGS/FigS13.eps','-depsc','-r300');

    end

    % here, instead of global/local, p-value based. For left-over scores a
    % different p-value score

    posShift = calc_pos_dif(oS,sortedIdsAll,synthStr{idxRun},lenThry, N);

    [passLocalThresh] = filter_scores_local(sortedVals,mpMaxLenBased,minOverlap,minLen,N);
    [passGlobalThresh] = filter_scores_global(foS,mpMaxLenBased,overlaplen,minLen,N);
    
    goodPos = passLocalThresh.*passGlobalThresh;
    import Plot.position_evaluation_plot;
    [fp] = position_evaluation_plot(passLocalThresh,passGlobalThresh,sortedVals,foS,posShift,N)

    % bad pos which are included in the model
    posShift(find(fp))

%     idxThreshMP = arrayfun(@(x) find(x<=minLen,1,'first'),minOverlap);
%     coefDiffsMP = sortedVals-mpMaxLenBased(idxThreshMP); % has to be already calculated
%     barsPassThreshMP = find(coefDiffsMP > 0);
%     
%     thresCC = mpMaxLenBased(idxThreshMP);
%     [sortedValsGood,sortedIdsGood,fullScoresGood] = filter_good_scores(oS,mpMaxLenBased,minLen,thresCC);

% thresCC = 0.9;
% import Core.filter_good_scores
% [sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMax,minLen,thresCC);

if plotfigs

% create a plot x:localscorethresh, y:numberof false positives
% figure

pxShiftLim = 20;
thresCCvals = 0.5:0.01:0.9;
fdr = zeros(1,length(thresCCvals));
nE= zeros(1,length(thresCCvals));
    import Core.filter_good_scores

for ii=1:length(thresCCvals)
    [sortedValsGood,sortedIdsGood,fullScoresGood,posShiftSel] = filter_good_scores(oS,mpMaxLenBased,minLen,thresCCvals(ii),posShift);
    fdr(ii) = sum(posShiftSel>=pxShiftLim)/length(posShiftSel);
    nE(ii) = length(posShiftSel);
end
figure;
hold on
tiledlayout(2,1);
nexttile
plot(thresCCvals,fdr);
title('fdr')

nexttile
plot(thresCCvals,nE);
title('number of elts passing threshold')
end

%%
if plotfigs

import Core.filter_good_scores
[sortedValsGood,sortedIdsGood,fullScoresGood] = filter_good_scores(oS,mpMaxLenBased,minLen,thresCC);


%   passthreshC{idxRun}  these passed thresh
figure,tiledlayout(2,1)
nexttile
plot(posShift,'x');hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
legend({'Distance between estimated GT and pair overlap','Threshold for good scores'},'Location','southoutside')
nexttile
plot(posShift(1:100));hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
legend({'Distance between estimated GT and pair overlap','Threshold for good scores'},'Location','southoutside')
print('FIGS/FigS5.eps','-depsc','-r300');



end


if plotfigs

    %
idxRun = 10;
import Core.comparisonsres_to_struct; 
[rezStruct2] = comparisonsres_to_struct(bG{idxRun}, comparisonStruct,theoryStruct{idxRun}{1}.length);
% 
%%
oS = resRun{idxRun}.oS;

idFold = 1;
idxPair = 5;

sets.minOverlap = 300;
 
% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, 0.08,0.01); %localCCStruct

[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
% 
% % lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
% lenThry = 10000;
% import Plot.pair_evaluation_plot;
% [f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,comparisonStructC{idxRun}([curSink curSource]),lenThry);



lenThry = theoryStruct{idxRun}{1}.length;
import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
[f] = pair_evaluation_with_ground_truth_plot(bG{idxRun}', oS,curSink,curSource,rezStruct2,lenThry);

% [f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,comparisonStructC{idxRun},lenThry);
print('FIGS/Fig1b.eps','-depsc','-r300');

      
sortedIdsGood = sortedIds(~isnan(sortedVals));
sortedValsGood = sortedVals( ~isnan(sortedVals));
import Validation.calc_pos_dif;
posShift = calc_pos_dif(oS,sortedIdsGood,comparisonStruct,lenThry, length(sortedIdsGood));

import Plot.position_evaluation_plot_stouffer;
[fp] = position_evaluation_plot_stouffer(sortedValsGood,localScore,posShift)

print('FIGS/FigS14.eps','-depsc','-r300');


end
% % pair_evaluation_with_ground_truth_plot
%     % We'll get a threshold for scores to get a positional difference plot.
% %     
% % idxThreshMP = arrayfun(@(x) find(x<=minLen,1,'first'),minOverlap);
% % coefDiffsMP = sortedVals-mpMaxLenBasedC{idxRun}(idxThreshMP); % has to be already calculated
% % barsPassThreshMP = find(coefDiffsMP > 0);
% % 
% % thresCC = mpMaxLenBasedC{idxRun}(idxThreshMP);
% % import Core.filter_good_scores
% % [sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBasedC{idxRun},minLen,thresCC);
% 
% length(sortedVals)
% 
% %%
% if plotfigs
% idFold = 1;
% idxPair = 1;
%  
% 
% 
% % allCoefs = cellfun(@(x) x.maxcoef(1),compStr);
% % allLengths =300;
% % 
% % coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% % 
% 
% 
% %
% 
% [curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
% 
% % lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
% lenThry = 10000;
% import Plot.pair_evaluation_plot;
% [f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,comparisonStructC{idxRun}([curSink curSource]),lenThry);
% print('FIGS/Fig1.eps','-depsc','-r300');
% 
% 

end



