
% get best indices
[a,idx] = max(cellfun(@(x) size(x,1),outConsensus));
% for i=1:size(outConsensus{idx},1)
%     plot((outConsensus{idx}(i,:)-nanmean(outConsensus{idx}(i,:)))/nanstd(outConsensus{idx}(i,:))+5*i,'black');
% end

% sort based on starting position
[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% 
[a,b] = max(cellfun(@(x) length(x.members),barsetGen.barcodeIslands));
% 
%
f = figure
g =tiledlayout(4,2,'TileSpacing','compact')
% nexttile
% hold on
% for i=1:length(barcodeGen)
%     plot(zscore(barcodeGen{i}.rawBarcode)+5*i,'black');
% end
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% title('(A) Best overlap')


idxPair = 84;
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));

lenThry =size(outConsensus{idx},2);
import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,[],lenThry,g);


nexttile([2 1])
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

highlight(h, indicesCell, 'NodeColor', 'r');
% 
% indicesCell1 = arrayfun(@(x) find(allNodes == x),  barsetGen.barcodeIslands{b}.edges(:,2), 'UniformOutput', true);
% indicesCell2 = arrayfun(@(x) find(allNodes == x),  barsetGen.barcodeIslands{b}.edges(:,3), 'UniformOutput', true);
% 
% 
% % highlightEdges = edges(GtempOrig, indicesCell', indicesCell');
% for i=1:length(indicesCell1)
%     highlight(h, 'Edges', [indicesCell1(i) indicesCell2(i)], 'EdgeColor', 'r');
% end
title('(B) Graph of barcode islands ','Interpreter','latex')

nexttile([2 1])

hold on

title('(C) Barcode island','Interpreter','latex')

consensusToPlot = outConsensus{idx}(idxv,:);
% B = nan(size(consensusToPlot,1) + size(consensusToPlot,1)-1, size(consensusToPlot, 2));

% Populate the new matrix with data from the original matrix
% B(1:2:end, :) = consensusToPlot;

imagesc(consensusToPlot);colormap(gray)
xlim([1 size(consensusToPlot,2)])
ylim([1 size(consensusToPlot,1)])

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

% 
nexttile([1 2])
% 
hold on

title('(D) Consensus barcode','Interpreter','latex')


imagesc(consensus{idx})
set(gca,'ytick',[])

xlim([1 size(consensusToPlot,2)])

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
xlabel('Position (kb)')

% 
print('FIGS/FigC2.eps','-depsc','-r300');

% next

%%


% comparisonStructCAll = horzcat( comparisonStructC{:} );
% thrSAll = horzcat( thrS{:} );

% import Core.comparisonsres_to_struct; 
% [rezStruct2] = comparisonsres_to_struct(barcodeGen, comparisonStructCAll,thrS{1}.length);
% 
% 
% import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
% f = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,comparisonStructCAll,lenThry);
% 
