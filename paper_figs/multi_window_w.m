%% OVERLAPS
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sF = 0.9:0.01:1.1; % length re-scaling factor
minOverlap = 200; % minimum overlap

plotfigs=0;

% can one run pre-processing procedure to reduce number of overlaps we need
% to calculate?


resRun = cell(1,length(bG));
idxRun = 1;
% for idxRun = 1:length(bG)
%     barcodeGen = bG{idxRun};
    lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen2);
    barcodeGen2 = barcodeGen2(lengths>minOverlap);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(barcodeGen2,sF, minOverlap, timestamp);
    toc
    resRun{idxRun}.oS = oS;
% end

save([num2str(dataSetIdx),'_',num2str(minOverlap),'_',timestamp,'resRun.mat'] ,'resRun','barcodeGen2' )

%% Generate results for each different window width w:
minO = [100:20:380];

for idxO = 1:length(minO)
    minOverlap = minO(idxO);

    resRun = cell(1,length(bG));
    
%     for idxRun = 1:length(bG)
    %     barcodeGen = bG{idxRun};
        barcodeGen3 = barcodeGen2(lengths>minOverlap);

        lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen2);
        tic
        % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
        [oS] = calc_overlap_mp(barcodeGen3,sF, minOverlap, timestamp);
        toc
        resRun{idxRun}.oS = oS;
%     end
    
    save([num2str(dataSetIdx),'_',num2str(minOverlap),'_',timestamp,'resRun.mat'] ,'resRun','barcodeGen3' )
end


%%


%%
oS = resRun{1}.oS;

alphaNu = 0.08;
pthresh = [0.001:0.001:0.01];
nB = 200;
sets.minOverlap = resRun{1}.oS(1,2).h;
sets.scDiffSetting = 0.04;
sets.pxDifSetting = 50;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

% ixtest = 5; % load the synthetic dataset with 3 islands

%%
pval =  0.01;
% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pval); %localCCStruct

sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

% sortedValsGood(4190:end) = nan;
% sortedIdsGood(1)
import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen3,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],1);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen2,[],1);


%%



% f=figure,g=tiledlayout(6,2,'TileSpacing','tight')
% pval =  pthresh(1);
% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pval); %localCCStruct
sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen3,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],0);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen3,[],0);

%%
f=figure,g=tiledlayout(4,2,'TileSpacing','loose')

idxPair = 2;
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
import Plot.pair_evaluation_with_ground_truth_simple_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_with_ground_truth_simple_plot(barcodeGen3, oS,curSink,curSource,[],[],g);


nexttile([1 2])

% plot(sortedVals,'o','Color',[0 0.4470 0.7410])
hold on
histogram(sortedVals,'EdgeAlpha',0)
histogram(sortedValsBad,'EdgeAlpha',0)
grid on
grid minor
lgd = legend({'$s_{combined, good}$ ','$s_{combined, bad}$ '},'Interpreter','latex')
lgd.Location = 'eastoutside';
set(gca, 'YScale', 'log')
title(['(C) Stouffer scores,pval=', num2str(pval)],'Interpreter','latex')




nexttile([1 2])
idx = 1;
% figure
totalDegree = indegree(GtempOrig) + outdegree(GtempOrig);
% Find nodes with total degree less than 2
nodesToRemove = find(totalDegree < 2);
G_removed = rmnode(GtempOrig, nodesToRemove);

h=plot(G_removed,'Layout','force','ArrowSize',5,'MarkerSize',1);
% 
allNodes = 1:numnodes(GtempOrig);
allNodes(nodesToRemove) = [];
indicesCell = arrayfun(@(x) find(allNodes == x),  barIslands{idx}, 'UniformOutput', true);

% highlight(h, indicesCell, 'NodeColor', 'r');
% 
title(['(C) Graph representation'], 'Interpreter','latex')

%

maxConsensus = max(cellfun(@(x) size(x,2),outConsensus));
% hold on

letters = {'D','E','F'}
for idx=1:length(outConsensus)
    nexttile([1 1])
hold on
title(['(' letters{idx} ') Block representation (', num2str(idx),')'], 'Interpreter','latex')

% idx = 1;

% sort based on starting position
[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% 
consensusToPlot1 = outConsensus{idx}(idxv,:);

imagesc(consensusToPlot1);colormap(gray)
% xlim([1 size(consensusToPlot1,2)])
ylim([1 size(consensusToPlot1,1)])

imagesc(outConsensus{idx}(idxv,:));colormap(gray)
axis off

set(gca,'xtick',[])

xlim([1 maxConsensus])
end
%

%     print('FIGS/Fig6.eps','-depsc','-r300');


