% Results for real data

load('3all_2023-11-24_12_38_55.mat')
load('3_2023-11-22_11_34_44resRun.mat')
load('3_individualdays.mat')
% final figure_6: Islands for dEC data:
%%

alphaNu = 0.08;
pthresh = [0.001:0.001:0.01];
nB = 200;
sets.minOverlap = 300;
sets.scDiffSetting = 0.04;
sets.pxDifSetting = 50;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

ixtest = 5; % load the synthetic dataset with 3 islands


nBvec=50:50:500;



f = figure;
tiledlayout(length(nBvec),length(pthresh));

ii=1;
for ii=1:length(nBvec);
for idxRun = 1:length(pthresh);
    nB = nBvec(ii);

%     oS = resRun{idxRun}.oS;
%     barcodeGen = bG{idxRun}';

    oSTemp = oS(1:nB,1:nB);
    bgTemp = barcodeGen(1:nB);
    
    % both on local and global
    [sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
        calculate_sorted_pvals(oSTemp,sets.minOverlap, 0.08,pthresh(idxRun)); %localCCStruct
    
    
    sortedValsGood = sortedVals(~isnan(sortedVals));
    sortedIdsGood  = sortedIds(~isnan(sortedVals));
    localScore = localScore(~isnan(sortedVals));

    
    import Core.barcode_island_output;
    [barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
        barcode_island_output(sortedValsGood,sortedIdsGood, oSTemp, bgTemp,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],0)
    %     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty
    
    % nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty
    
        import Core.create_graph_from_barsetgen;
     [Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,bgTemp,[],0);


    
    nexttile
    % figure
    totalDegree = indegree(GtempOrig) + outdegree(GtempOrig);
    % Find nodes with total degree less than 2
    nodesToRemove = find(totalDegree < 2);
    G_removed = rmnode(GtempOrig, nodesToRemove);
    
    h=plot(G_removed,'Layout','force','ArrowSize',5,'MarkerSize',1,'NodeLabel', {});
    % 
    allNodes = 1:numnodes(GtempOrig);
    allNodes(nodesToRemove) = [];
    %     indicesCell = arrayfun(@(x) find(allNodes == x),  barIslands{idx}, 'UniformOutput', true);
        
    %     highlight(h, indicesCell, 'NodeColor', 'r');
    if ii==1
        title(['SNR = ', num2str(idxRun)],'Interpreter','latex')
    end
    
    ax = gca;
    
    % Rotate y-axis labels vertically
    set(ax, 'YTickLabelRotation', 90);
    if idxRun==2
    ylabel(['#barcodes = ',num2str(nB)],'FontSize',6);
    end
    
    end
end

print('FIGS/FigTemp6_1.eps','-depsc','-r300');


%%



% maybe put into sep function
% oS = resRun{ixtest}.oS;
% barcodeGen = bG{ixtest}';

% % oS = oS(1:nB,1:nB);
% barcodeGen = barcodeGen(1:nB);

pval =  pthresh(end);
% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pval); %localCCStruct

sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],1);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen,[],1);



f=figure,g=tiledlayout(6,2,'TileSpacing','tight')
nexttile([1 2])
pval =  pthresh(1);
% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pval); %localCCStruct


% plot(sortedVals,'o','Color',[0 0.4470 0.7410])
hold on
histogram(sortedVals,'EdgeAlpha',0)
histogram(sortedValsBad,'EdgeAlpha',0)
grid on
grid minor
lgd = legend({'$s_{combined, good}$ ','$s_{combined, bad}$ '},'Interpreter','latex')
lgd.Location = 'eastoutside';
set(gca, 'YScale', 'log')
title(['(A) Stouffer scores,pval=', num2str(pval)],'Interpreter','latex')
nexttile([1 2])
pval =  pthresh(end);
% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pval); %localCCStruct


% plot(sortedVals,'o','Color',[0 0.4470 0.7410])
hold on
histogram(sortedVals,'EdgeAlpha',0)
histogram(sortedValsBad,'EdgeAlpha',0)
grid on
grid minor
set(gca, 'YScale', 'log')
lgd = legend({'$s_{combined, good}$ ','$s_{combined, bad}$ '},'Interpreter','latex')
lgd.Location = 'eastoutside';
title(['(B) Stouffer scores,pval=', num2str(pval)],'Interpreter','latex')




[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pval); %localCCStruct

sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],0);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen,[],0);



maxConsensus = max(cellfun(@(x) size(x,2),outConsensus));
% hold on

letters = {'C','D','E','F'}
for idx=1:4
    nexttile([1 2])
hold on
title(['(' letters{idx} ') Aligned barcode island (', num2str(idx),')'], 'Interpreter','latex')

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


    print('FIGS/Fig6test1.eps','-depsc','-r300');
