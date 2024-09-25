% load data. Here it changes to folder with external data. For
% reproducibility, the '.mat' file with the overlaps can be downloaded
cd('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/')

%
load('3all_2023-11-24_12_38_55.mat','oS','barcodeGen')
% load('3_2023-11-22_11_34_44resRun.mat')
% load('3_individualdays.mat')

% load('experiment_data_oS.mat');
% load('experiment_data_barcodeGen.mat');


%% Parameters for the run. Check that these to be the same as Table 2 in the paper


% current in paper: 0.075, 0.001
alphaNu1 = 0.085; % 0.08 default for the paper
alphaNu2 = 0.065; % 0.08 default for the paper

alphaN1 = 0.15; % 0.08 default for the paper
alphaN2 = 0.002; % 0.08 default for the paper

sets.minOverlap = 300;
sets.scDiffSetting = 0.05;
sets.pxDifSetting = 50;
showFig = 1; % whether to show output figure
pthresh =  0.05;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

%
% both on local and global. Only uses the oS structure as input +
% parameters from table 2
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu1, pthresh,[],alphaN1,alphaNu2,alphaN2); %localCCStruct


sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

% sortedValsGood(4190:end) = nan; % change alphaNu instead of this

% Main barcode islands script:

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],showFig,50);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty
accepted = cellfun(@(x) x.accepted, barsetGen.currentData); % all the

import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen,[],showFig);

%% Intro scripts:
% intro_rescaling.mlx 

wminC = 300;
NN = 5;

import Core.barcode_island_consensus;
aaBarcodeIslandBlockRep = cell(1,length(cGenAll));
clusterBarcodes =  cell(1,length(cGenAll));
matRep =  cell(1,length(cGenAll));
consensusNew =  cell(1,length(cGenAll));
barCon =  cell(1,length(cGenAll));
stats =  cell(1,length(cGenAll));
cIt =  cell(1,length(cGenAll));

parfor ix = 1:length(cGenAll)
%     ix
    [aaBarcodeIslandBlockRep{ix}, ~,~,clusterBarcodes{ix},aaBarcodeIslandBlockRepRaw{ix}] = barcode_island_consensus(barcodeGen,cGenAll, ix, wminC);

    % wminC
    [matRep{ix},consensusNew{ix},barCon{ix},stats{ix},cIt{ix}] = iterative_readjustment_consensus(aaBarcodeIslandBlockRep{ix}, barcodeGen, barIslands, cGenAll, ix, wminC, NN,alphaNu1);

end


outConsensus = cellfun(@(x) x{end},matRep,'un',false);

save('experiment_matisland.mat','matRep','consensusNew','aaBarcodeIslandBlockRep','aaBarcodeIslandBlockRepRaw','clusterBarcodes','barCon','stats','cIt','barsetGen', 'outConsensus', 'coverage', 'consensus', 'islandElts', 'islandPx','cGenAll','barcodeIslandsData', 'barStruct','barIslands');

%% Here also plot for FigS8:
% consensus_real_data.mlx
ix = 2;
 f = consensus_ci(aaBarcodeIslandBlockRep{ix},aaBarcodeIslandBlockRepRaw{ix},matRep{ix}{end})
print('FIGS/Fig6SS2.eps','-depsc','-r500');

ix = 1;
 f = consensus_ci(aaBarcodeIslandBlockRep{ix},aaBarcodeIslandBlockRepRaw{ix},matRep{ix}{end})
print('FIGS/Fig6SS1.eps','-depsc','-r500');

ix = 3;
 f = consensus_ci(aaBarcodeIslandBlockRep{ix},aaBarcodeIslandBlockRepRaw{ix},matRep{ix}{end})
print('FIGS/Fig6SS3.eps','-depsc','-r500');
%% Figure6 in the paper:

f=figure('Position', [10 10 714 800]); g=tiledlayout(6,2,'TileSpacing','tight','Padding','none');
% fig = gcf;
% fig.PaperSize=[21 29.7];
% fig.PaperPosition = [1 1 20 25];

idxPair = 69;
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
import Plot.pair_evaluation_with_ground_truth_simple_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,curSink,curSource,[],[],g);


nexttile([1 2]) % Stouffer scores, Good, bad, & scores that are put into islands

% plot(sortedVals,'o','Color',[0 0.4470 0.7410])
hold on
histogram(sortedVals,'EdgeAlpha',0.1,'BinLimits',[-5 10],'BinWidth',0.5)
histogram(sortedValsBad,'EdgeAlpha',0.1,'BinLimits',[-5 10],'BinWidth',0.5)

[N,edges] = histcounts(sortedValsGood(logical(accepted)),'Normalization','count','BinLimits',[-5 10],'BinWidth',0.5);
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N,'Color','black');

% histogram(,'EdgeAlpha',1,'BinLimits',[-5 10],'BinWidth',0.5, 'DisplayStyle', 'stairs')

grid on
grid minor
lgd = legend({'$s_{\rm Stouffer}$ ','$s_{\rm Removed}$ ','$s_{Final}$'},'Interpreter','latex')
lgd.Location = 'eastoutside';
set(gca, 'YScale', 'log')
title(['(B) Stouffer scores'],'Interpreter','latex')

totalDegree = indegree(GtempOrig) + outdegree(GtempOrig);
% Find nodes with total degree less than 2
nodesToRemove = find(totalDegree < 2);
allNodes = 1:numnodes(GtempOrig);
allNodes(nodesToRemove) = [];
G_removed = rmnode(GtempOrig, nodesToRemove);

lettersAll = mat2cell('A':'Z',1,ones(1,26));

letters = lettersAll(3:end);
nexttile([1 2]) % graphs
hold on
idx=1;
% title(['(' letters{idx} ') Graph representation (', num2str(idx),')'], 'Interpreter','latex')
title(['(' letters{idx} ') Graph representation '], 'Interpreter','latex')

for idx=1:length(outConsensus)
    indicesCell = arrayfun(@(x) logical(allNodes == x),  barIslands{idx}, 'UniformOutput', false);

    allNodesKeep = sum(cell2mat(indicesCell'));
    % remove all except these
    G_representation = rmnode(G_removed, find(~allNodesKeep));

    h=plot(G_representation,'Layout','force','ArrowSize',5,'MarkerSize',1);
% 
    h.NodeLabelMode = 'auto'; 
    % To remove all labels
    h.NodeLabel = {};
    h.XData = h.XData+10*idx;
    text(min(h.XData),min(h.YData)-1, [ '(' num2str(idx) ')'])
end
%
axis off

maxConsensus = max(cellfun(@(x) size(x,2),outConsensus));
% hold on

letters = lettersAll(4:end)
for idx=1:length(outConsensus)
    nexttile([1 2])
    hold on
    title(['(' letters{idx} ') Block representation (', num2str(idx),') zoom-in (center 2000 pixels)'], 'Interpreter','latex')
    
    % sort based on starting position
    [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
    % 
    consensusToPlot1 = outConsensus{idx}(idxv,:);
    
    % zoomin
    consensusToPlot1 = consensusToPlot1(:,[round(size(consensusToPlot1,2)/2-1000:size(consensusToPlot1,2)/2+1000-1)]);
    consensusToPlot1 = consensusToPlot1(find(sum(~isnan(consensusToPlot1'))),:);
    
    consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/5),1);
    consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
%     imagesc([consensusToPlot1;consensusEmptyBlock;consensusBlock]);colormap(gray)
    imagesc([consensusBlock;consensusEmptyBlock;consensusToPlot1;]);colormap(gray)

%     set(gca,'YDir','normal')
    ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
    text(size(consensusBlock,2)+10,size(consensusBlock,1)/2, 'Consensus','FontSize',7, 'FontName','Times',   'Color','black')

axis off

set(gca,'xtick',[])

xlim([1 2000])

end
%
nPixels = 1e5/500;% 500bp/px ;
x = [2000-nPixels-100 2000-100];
y = [20.6 20.6 ];
plot(x,35+y,'Linewidth',2,'Color','white')
% text(x(1),y(1)+0.1, '10 $\mu m$','FontSize',12, 'Color','white','Interpreter','latex')
text(x(1)+nPixels/10,y(1)+52, '100 kb','FontSize',14, 'FontName','Times',   'Color','white')
set(gcf, 'Color', 'w')

f.Position(3) = 714;

print('FIGS/Fig6.eps','-depsc','-r500');


    %%

    for idx=1:length(outConsensus)
        figure;
        ax = gca;

        % sort based on starting position
        [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
        %
        consensusToPlot1 = outConsensus{idx}(idxv,:);

        consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/5),1);
        consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
        %     imagesc([consensusToPlot1;consensusEmptyBlock;consensusBlock]);colormap(gray)
        imagesc([consensusToPlot1; consensusEmptyBlock;consensusBlock]);colormap(gray)

        %     set(gca,'YDir','normal')
        ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
        %     axis off
        set(gcf, 'Color', 'w')
        hold on

        nPixels = 1e5/500;% 500bp/px ;
        x = [size(consensusToPlot1,2)-460 size(consensusToPlot1,2)-460+nPixels];
        y = [2 2 ];
        plot(x,y,'Linewidth',2,'Color','white')
%         text(size(consensusBlock,2)+10,size(consensusToPlot1,1)+80, 'Consensus','FontSize',7, 'FontName','Times',   'Color','black')

        ax.YAxis.Visible = 'off';
        xlabel('Position (px)','Interpreter','latex')
        print(['FIGS/Fig6S',num2str(idx) , '.eps'],'-depsc','-r500');

    end
