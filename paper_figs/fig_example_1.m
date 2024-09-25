% Example Figure 1: 

% PAPER DATA IN /export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/

% Schematics of de-novo assembly problem using DNA barcodes.
rng('default')

snr = 5; 
import Core.load_synth_data;
[bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(30,2.72,850,100,snr,0.05,1,2000,150);


% we have pre-generated, so just load pre-generated data
% load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/synth_2023-11-28_09_43_54bml.mat');
% load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/synth_2023-11-28_09_43_54resRun.mat');

sets.minOverlap = 300;
timestamp ='';
idxRun = 1;
oS = synthStr2{idxRun};%resRun{idxRun}.oS;
barcodeGen = bG{idxRun}';

% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS,nan, 0.08,0.01); %localCCStruct

sortedVals = sortedVals(1:length(localScore));
sortedVals(100:end) = nan; % keep less scores to show two islands
sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.02;
sets.pxDifSetting = 20; %20

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [])
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

import Core.create_graph_from_barsetgen;
[Gtemp] = create_graph_from_barsetgen(barsetGen,barcodeGen);



%
f = figure
tiledlayout(2,2,'TileSpacing','compact')
nexttile
hold on
for i=1:length(barcodeGen)
    plot(zscore(barcodeGen{i}.rawBarcode)+5*i,'black');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(A) List of barcodes')

nexttile

plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1);
title('(B) Overlap graph ')

nexttile([1 2])
hold on

title('(C) Largest barcode island ', 'Interpreter','latex')

[a,idx] = max(cellfun(@(x) size(x,1),outConsensus));
for i=1:size(outConsensus{idx},1)
    plot((outConsensus{idx}(i,:)-nanmean(outConsensus{idx}(i,:)))/nanstd(outConsensus{idx}(i,:))+5*i,'black');
end
set(gca,'ytick',[])
xlabel('Position (kb)')

xticksOriginal = get(gca, 'XTick');
xtickLabelsOriginal = get(gca, 'XTickLabel');

scaleFactor = 0.5;
% Scale the tick locations
xticksScaled = xticksOriginal * scaleFactor;


% Optionally, update the tick labels (if you want them to reflect the scaled values)
xtickLabelsScaled = cellstr(num2str(xticksScaled'));
set(gca, 'XTickLabel', xtickLabelsScaled);

print('FIGS/Fig1.eps','-depsc','-r300');

% nexttile
% hold on
% 
% title('(D) Barcode island 2 ')


%% Example to table (for fig2B):


%
f = figure
tiledlayout(4,2,'TileSpacing','compact','Padding','none')
nexttile([2 1])
hold on
for i=1:length(barcodeGen)
    plot(zscore(barcodeGen{i}.rawBarcode(barcodeGen{i}.rawBitmask))+5*i,'black');
end
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(A) List of barcodes', 'Interpreter','latex')

nexttile([2 1])

plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1);
title('(B) Graph of barcode islands ', 'Interpreter','latex')

% nexttile([1 2])
% hold on
% 
% title('(C) Barcode island ')
% 
% [a,idx] = max(cellfun(@(x) size(x,1),outConsensus));
% for i=1:size(outConsensus{idx},1)
%     plot((outConsensus{idx}(i,:)-nanmean(outConsensus{idx}(i,:)))/nanstd(outConsensus{idx}(i,:))+5*i,'black');
% end
% set(gca,'ytick',[])
% xlabel('Position (kb)')
% 
% xticksOriginal = get(gca, 'XTick');
% xtickLabelsOriginal = get(gca, 'XTickLabel');
% 
% scaleFactor = 0.5;
% % Scale the tick locations
% xticksScaled = xticksOriginal * scaleFactor;
% 
% 
% % Optionally, update the tick labels (if you want them to reflect the scaled values)
% xtickLabelsScaled = cellstr(num2str(xticksScaled'));
% set(gca, 'XTickLabel', xtickLabelsScaled);



nexttile([1 2])

hold on
title('(C) Aligned barcode island', 'Interpreter','latex')

% get best indices
[a,idx] = max(cellfun(@(x) size(x,1),outConsensus));

% sort based on starting position
[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% 
consensusToPlot = outConsensus{idx}(idxv,:);

imagesc(consensusToPlot);colormap(gray)
xlim([1 size(consensusToPlot,2)])
ylim([1 size(consensusToPlot,1)])

imagesc(outConsensus{idx}(idxv,:));colormap(gray)
set(gca,'xtick',[])
nexttile([1 2])

hold on
title('(D) Consensus barcode', 'Interpreter','latex')

imagesc(consensus{idx})
set(gca,'ytick',[])
xlabel('Position (kb)')

xticksOriginal = get(gca, 'XTick');
xtickLabelsOriginal = get(gca, 'XTickLabel');% nexttile([1 2])


% xticksOriginal = get(gca, 'XTick');
% xtickLabelsOriginal = get(gca, 'XTickLabel');

scaleFactor = 0.5;
% Scale the tick locations
xticksScaled = xticksOriginal * scaleFactor;
% 
% 
% Optionally, update the tick labels (if you want them to reflect the scaled values)
xtickLabelsScaled = cellstr(num2str(xticksScaled'));
set(gca, 'XTickLabel', xtickLabelsScaled);
% 
xlim([1 size(consensusToPlot,2)])
% ylim([1 size(consensusToPlot,1)])

print('FIGS/Fig1.eps','-depsc','-r300');
