%% OVERLAPS
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sF = 0.9:0.01:1.1; % length re-scaling factor
minOverlap = 300; % minimum overlap
minOverlapContig = 600;


plotfigs=0;

% can one run pre-processing procedure to reduce number of overlaps we need
% to calculate?



resRun = cell(1,length(bG));

for idxRun = 1:length(bG)
    barcodeGen = bG{idxRun};
    lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(barcodeGen,sF, minOverlap, timestamp);
    toc
    resRun{idxRun}.oS = oS;
end

save([num2str(dataSetIdx),'_',timestamp,'resRun.mat'] ,'resRun','bG' )


% this could be calculated via _individual_days
try
nmPerPx = kymoStructs{idxRun}{1}.nmpxidFold;
nmbp = kymoStructs{idxRun}{1}.nmBpidFold;
catch
    warning('params not provided')
    nmPerPx = 110;
    nmbp = 0.25;
end

%%--- PVALS
numWorkers = 30;
fastaFile = fastas((theoryIdxs{dataSetIdx}));
minLen = [150:50:3000]; % min to max % could take more points for more accurate..
% can be calculated before for all nmbp
[MP,mpMaxLenBasedC{idxRun},theoryStructRev,MPI,thrS{idxRun}] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,1,numWorkers,minLen, []);% no lenseq if based on fastafile
%
localCCStruct.minLen = minLen;
localCCStruct.maxBasedCC = mpMaxLenBasedC{idxRun};



%     barcodeGen = barcodeGen

% 
% % 
% % idxs = 1:length(bgAll);
% % barcodeGen = bgAll(barsPassThresh);idxs = barsPassThresh;
% 
% bgS = barcodeGen;
% 
% % [oS] = calc_scores({bgAll},300,0.95:0.025:1.05);
% 
% tic
% % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
% [oS] = calc_overlap_mp(bgS,sF, minOverlap,timestamp);
% toc
% 
%% bargrouping_funs

idxRun = 9;
oS = resRun{idxRun}.oS;
barcodeGen = bG{idxRun};

% both on local and global
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
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

    import Core.create_graph_from_barsetgen;
 [Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen);



% validation against theory
import Validation.barislands_vs_theory_validation

barislands_vs_theory_validation

import Validation.barislands_vs_theory_validation_individual

barislands_vs_theory_validation_individual
% 
% import Core.plot_match_pcc;
% [sortedVals, sortedIdsAll,foS] = sorted_scores(oS);
% idxThreshMP = arrayfun(@(x) find(x<=minLen,1,'first'),minOverlap);
% coefDiffsMP = sortedVals-mpMaxLenBasedC{idxRun}{1}(idxThreshMP); % has to be already calculated
% barsPassThreshMP = find(coefDiffsMP > 0);
% 
% thresCC = mpMaxLenBasedC{idxRun}{1}(idxThreshMP);
% import Core.filter_good_scores
% [sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBasedC{idxRun}{1},minLen,thresCC);

length(sortedVals(~isnan(sortedVals)))

%%
if plotfigs

    %
import Core.comparisonsres_to_struct; 
[rezStruct2] = comparisonsres_to_struct(barcodeGen, comparisonStructC{idxRun},thrS{idxRun}.length);
% 
%%
idFold = 1;
idxPair = 5;
 

[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
% 
% % lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
% lenThry = 10000;
% import Plot.pair_evaluation_plot;
% [f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,comparisonStructC{idxRun}([curSink curSource]),lenThry);



lenThry = thrS{idxRun}.length;
import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
[f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,rezStruct2,lenThry);

% [f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,comparisonStructC{idxRun},lenThry);
print('FIGS/Fig1b.eps','-depsc','-r300');

      
import Validation.calc_pos_dif;
posShift = calc_pos_dif(oS,sortedIdsGood,comparisonStructC{idxRun},lenThry, length(sortedIdsGood));

import Plot.position_evaluation_plot_stouffer;
[fp] = position_evaluation_plot_stouffer(sortedValsGood,localScore,posShift)

%         lenThry = theoryStruct{idxRun}{1}.length;
%         [f] = pair_evaluation_with_ground_truth_plot(barcodeGen', oS,curSink,curSource,synthStr2{idxRun},lenThry);
  

%%
N = 1000;
posShift = zeros(1,N);
for idxPair=1:N
    
    [curSink,curSource] = ind2sub(size(oS),sortedIdsAll(idxPair));
    pos1 = comparisonStructC{idxRun}{curSink}.pos(1);
    pos2 = comparisonStructC{idxRun}{curSource}.pos(1);
    if comparisonStructC{idxRun}{curSource}.or(1)==1
        posDifTheory = pos1-pos2;
    else
        posDifTheory = -(pos1-pos2+comparisonStructC{idxRun}{curSink}.lengthMatch-comparisonStructC{idxRun}{curSource}.lengthMatch);
    end
    posDifExp = oS(curSink,curSource).pB - oS(curSink,curSource).pA;


    
    posShift(idxPair) = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);
end

%   passthreshC{idxRun}  these passed thresh
figure,tiledlayout(2,1)
nexttile
plot(posShift);hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})
nexttile
plot(posShift(1:100));hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})
print('FIGS/FigS3.eps','-depsc','-r300');

import Validation.calc_pos_dif;
import Plot.pair_evaluation_with_ground_truth_plot;
import Validation.calc_pos_dif;
import Nullmodel.filter_scores_local; % all scores up to a value, since  minOverlap is fixed   
import Nullmodel.filter_scores_global; % all scores up to a value, since  minOverlap is fixed   

posShift = calc_pos_dif(oS,sortedIdsGood,comparisonStructC{idxRun},lenThry, N);

% [passLocalThresh] = filter_scores_local(sortedValsGood,mpMaxLenBasedAll{idxRun},minOverlap,minLen,N);
% [passGlobalThresh] = filter_scores_global(foS,mpMaxLenBasedAll{idxRun},overlaplen,minLen,N);

goodPos = passLocalThresh.*passGlobalThresh;
import Plot.position_evaluation_plot;
[fp] = position_evaluation_plot(passLocalThresh,passGlobalThresh,sortedValsGood,foS,posShift,N)

end
%%
% temp
[sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBasedC{idxRun}{1},minLen,thresCC);
ix=4;
sortedIds = sortedIds(1:ix);
sortedVals = sortedVals(1:ix);
[curSink,curSource] = ind2sub(size(oS),sortedIds(1));

%
% plot the graph
import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(size(oS),sortedIds,length(sortedIds),barcodeGen,zeros(1,length(pscores)),100,oS);
print('FIGS/FigS4.eps','-depsc','-r300');

%
% [mpMax,finalgraph,Ggraphs,sortedVals,sortedIds,pscores] = bargrouping_funs(oS,bG,0.8)
% toc
import Core.bar_island_out;
[outConsensus, coverage, consensus,islandElts, islandPx,cGen,badData,dataInfo,barcodeIslandsData,barcodeIslands,barStruct] = bar_island_out(sortedVals,sortedIds, oS,barcodeGen,timestamp);
print('FIGS/Fig2.eps','-depsc','-r300');

%%
import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmPerPx);
barcodeGenT{1}.rawBitmask = logical(barcodeGenT{1}.rawBitmask);

sets.comparisonMethod = 'mass_pcc';

%  add theory as last barcode
cGen =[cGen barcodeGenT];

sF = 0.9:0.01:1.1;
% [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(cGen,theoryStruct,sF,sets);
% [compStr,~,calcLengths] = compare_to_t_mp(cGen,theoryStruct,sF,300); % get oS as output

% calculate MP scores 
[oScons] = calc_overlap_mp(cGen,sF, minOverlapContig,timestamp);
 
import Core.plot_match_simple;
[f] = plot_match_simple(cGen, oScons, 8,10);
print('FIGS/Fig4.eps','-depsc','-r300');


% [f] = plot_match_simple(bgS, oS,curSink,curSource);


% 
% allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
% allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);
% 
% idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
% coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% 
% barsPassThresh = find(coefDiffs > 0.05);
% % barsPassThresh = find(ones(1,length(comparisonStruct)));

super_quick_plot(910,bgAll,comparisonStruct,theoryStruct)


compare_consensuses()


end
%% Now validation for full


%%

import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(bgS, oS,curSink,curSource);

bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct(idxs([curSink curSource])), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6


% plot on theory

% %% check position shift (ignoring orientation for now)
% N = 1000;
% posShift = zeros(1,N);
% for idxPair=1:N
%     
%     [curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
%     pos1 = comparisonStruct{idxs([curSink])}.pos(1);
%     pos2 = comparisonStruct{idxs([curSource])}.pos(1);
%     if comparisonStruct{idxs([curSource])}.or(1)==1
%         posDifTheory = pos1-pos2;
%     else
%         posDifTheory = -(pos1-pos2+comparisonStruct{idxs([curSink])}.lengthMatch-comparisonStruct{idxs([curSource])}.lengthMatch);
%     end
%     posDifExp = oS(curSink,curSource).pB - oS(curSink,curSource).pA;
% 
% 
%     
%     posShift(idxPair) = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);
% end
% 
% figure,plot(posShift);hold on
% plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
% legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})


%% bargrouping_funs
thresCC = mpMaxLenBased(idxThreshMP);
import Core.filter_good_scores
[sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBased,minLen,thresCC);

%%
idFold = 1;
idxPair = 7;
 


% allCoefs = cellfun(@(x) x.maxcoef(1),compStr);
% allLengths =300;
% 
% coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% 


%

[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));

% lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
lenThry = 10000;
import Plot.pair_evaluation_plot;
[f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,comparisonStruct([curSink curSource]),lenThry);

super_quick_plot(66,bgAll,comparisonStruct,theoryStruct)
