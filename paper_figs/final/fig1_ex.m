% Example Figure 1: 

% originally fig_example_1
%
% PAPER DATA IN /export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/

% Schematics of de-novo assembly problem using DNA barcodes.
rng('default')

snr = 5; 
import Core.load_synth_data;
[bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(30,2.72,850,100,snr,0.05,1,2000,150,0,[]);
% [bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(50,2.72,850,100,snr,0.05,1,2000,150);


% we have pre-generated, so just load pre-generated data
% load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/synth_2023-11-28_09_43_54bml.mat');
% load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/synth_2023-11-28_09_43_54resRun.mat');

sets.minOverlap = 300;
timestamp ='';
idxRun = 1;
oS = synthStr2{idxRun};%resRun{idxRun}.oS;
barcodeGen = bG{idxRun}';

% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver,pvalCombined,  sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,nan, 0.08,0.01); %localCCStruct

sortedVals = sortedVals(1:length(localScore));
sortedVals(100:end) = nan; % keep less scores to show two islands
% sortedVals(2000:end) = nan; % keep less scores to show two islands
sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.02;
sets.pxDifSetting = 20; %20

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] = ...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, []);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

import Core.create_graph_from_barsetgen;
[Gtemp] = create_graph_from_barsetgen(barsetGen,barcodeGen);

% get best indices
consSize = cellfun(@(x) size(x,1),outConsensus);

[sV,sI] = sort(consSize,'descend');
% [a,idx] = max(consSize);


%% Example to table (for fig2B):


%
f=figure('Position', [10 10 714 600]);
tiledlayout(65,2,'TileSpacing','tight','Padding','none')
nexttile([26 1])
hold on
for i=1:length(barcodeGen)
    plot(zscore(barcodeGen{i}.rawBarcode(barcodeGen{i}.rawBitmask))+5*i,'black');
    if i<3 || i >length(barcodeGen)-2
    text(-40,5*i,num2str(i),'FontSize',5)
    else
    text(-40,5*i,' :','FontSize',5)
    end
end
set(gca,'xtick',[])
set(gca,'ytick',[])
title('(A) List of barcodes', 'Interpreter','latex')

nexttile([26 1])

plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1);
title('(B) Graph representation', 'Interpreter','latex')

text(0,4,'(1)')
text(0,8,'(2)')
text(3.2,8,'(3)')
text(4.8,2.3,'(4)')


nexttile([15 2])

hold on
title('(C) Block representation of barcode island (1)', 'Interpreter','latex')

idx = sI(1);

% sort based on starting position
[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% 
consensusToPlot = outConsensus{idx}(idxv,:);

imagesc(consensusToPlot);colormap(gray)
xlim([1 size(consensusToPlot,2)])
% ylim([0 size(consensusToPlot,1)+1])

imagesc(outConsensus{idx}(idxv,:));colormap(gray)
set(gca,'xtick',[])

% title('(C) Aligned bargroup island (1)', 'Interpreter','latex')

import Plot.block_plot;

consensusToPlot1 = block_plot(sI,1, outConsensus);

% 
% imagesc(consensusToPlot1);colormap(gray)
% xlim([1 size(consensusToPlot1,2)])
% ylim([0.5 size(consensusToPlot1,1)+0.5])
% 
% imagesc(outConsensus{idx}(idxv,:));colormap(gray)
% set(gca,'xtick',[])
% axis off


% 
nexttile([7 2])

hold on
title('(D) Consensus barcode of barcode island (1)', 'Interpreter','latex')
idx = sI(1);

imagesc(consensus{idx})
set(gca,'ytick',[])
xlim([1 size(consensusToPlot1,2)])
axis off
% 
% 

nPixels = 1e5/500;% 500bp/px ;
x = [size(consensusToPlot1,2)-nPixels-100 size(consensusToPlot1,2)-100];
y = [0.6 0.6 ];
plot(x,y,'Linewidth',2,'Color','white')
% text(x(1),y(1)+0.1, '10 $\mu m$','FontSize',12, 'Color','white','Interpreter','latex')
text(x(1)+nPixels/5,y(1)+0.4, '100 kb','FontSize',14, 'FontName','Times',   'Color','white')
set(gcf, 'Color', 'w')

nexttile([10 2])
hold on
title('(E) Block representation of barcode island (2)', 'Interpreter','latex')

import Plot.block_plot;

consensusToPlot1 = block_plot(sI,2,outConsensus);

% % get best indices
% % [a,idx] = max(cellfun(@(x) size(x,1),outConsensus));
% idx = sI(2);
% 
% % sort based on starting position
% [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% % 
% consensusToPlot = outConsensus{idx}(idxv,:);
% 
% imagesc(consensusToPlot);colormap(gray)
% xlim([1 size(consensusToPlot1,2)])
% ylim([0.5 size(consensusToPlot,1)+0.5])
% 
% imagesc(outConsensus{idx}(idxv,:));colormap(gray)
% axis off
% set(gca,'xtick',[])

nexttile([7 2])

hold on
title('(F) Consensus barcode of barcode island (2)', 'Interpreter','latex')
idx = sI(2);

imagesc(consensus{idx})
set(gca,'ytick',[])
xlim([1 size(consensusToPlot1,2)])
axis off
% 
% xxlabel('Position (kb)')
% ticksOriginal = get(gca, 'XTick');
% xtickLabelsOriginal = get(gca, 'XTickLabel');% nexttile([1 2])

% scaleFactor = 0.5;
% % Scale the tick locations
% xticksScaled = xticksOriginal * scaleFactor;
% 
% xtickLabelsScaled = cellstr(num2str(xticksScaled'));
% set(gca, 'XTickLabel', xtickLabelsScaled);
% 

nPixels = 1e5/500;% 500bp/px ;
x = [size(consensusToPlot1,2)-nPixels-100 size(consensusToPlot1,2)-100];
y = [0.6 0.6 ];
plot(x,y,'Linewidth',2,'Color','white')
% text(x(1),y(1)+0.1, '10 $\mu m$','FontSize',12, 'Color','white','Interpreter','latex')
text(x(1)+nPixels/5,y(1)+0.4, '100 kb','FontSize',14, 'FontName','Times',   'Color','white')
set(gcf, 'Color', 'w')

print('FIGS/Fig1.eps','-depsc','-r500');

%%
extra = 0;
if extra
%%


N = length(sortedValsGood);

import Validation.calc_distance_scaled; % need to properly scale things when calculating distances for validation
posShift = calc_distance_scaled(oS,sortedIdsGood(1:N), synthStr{1},theoryStruct{1}{1}.length, N);

N = min(1000,length(sortedIdsBad));


import Validation.calc_distance_scaled;
posShiftBad = calc_distance_scaled(oS,sortedIdsBad(1:N), synthStr{1},theoryStruct{1}{1}.length, N);

N = length(sortedValsGood);

import Validation.calc_pos_dif;
posShift = calc_pos_dif(oS,sortedIdsGood(1:N), synthStr{1},theoryStruct{1}{1}.length, N);


N = min(1000,length(sortedIdsBad));

import Validation.calc_pos_dif;
posShiftBad = calc_pos_dif(oS,sortedIdsBad(1:N), synthStr{1},theoryStruct{1}{1}.length, N);

% % 
% import Plot.position_evaluation_stouffer;
% [f,g] = position_evaluation_stouffer(sortedValsGood,sortedValsBad, localScore,posShift,posShiftBad)
% 



% this function updates data vector. Still need to get barcodes for this
import Core.update_sf_barset;
barIslands =barIslands(1);
ix = 5;
 
pDif = cell(1,length(barIslands));
sorti =  cell(1,length(barIslands));
val =  cell(1,length(barIslands));
for idx = 1%:length(barIslands)

    
    bars = barcodeGen(barIslands{idx}); 
    barx = barIslands{idx};
    
    synCur = synthStr{1}(barx);


    bbS = synCur{ix}.bestBarStretch;  % best bar stretch (bbS)
    bbO = (synCur{ix}.or(1)~=barcodeIslandsData{idx}(ix,3))+1;


    [bbSAll] = update_sf_barset(barcodeIslandsData{idx}, bbS/barcodeIslandsData{idx}(ix,4), bbO);

    % need to also update positions for all the synCur barcodes based on
    import Core.synth_to_table;
    [tableS] = synth_to_table(synCur);

    % update synth table. Want to have the same
    [tableSUpd] = update_sf_barset(tableS, tableS(ix,4)/bbS, 1);


    bbSAll(:,1:2) = bbSAll(:,1:2)-bbSAll(ix,1)+tableSUpd(ix,1); %
   
    bestPosFound = bbSAll(:,1);
    
    posTrue = tableSUpd(:,1);
    %arrayfun(@(x) synthStr{1}{x}.pos, barx);
    
%     pDif{idx} = min([sqrt((posTrue'-bestPosFound).^2) sqrt((posTrue'+theoryStruct{1}{1}.length-bestPosFound).^2)...
%         sqrt((posTrue'-theoryStruct{1}{1}.length-bestPosFound).^2)]');


    pDif{idx} = posTrue-bestPosFound;

     pDif{idx}( pDif{idx}<-theoryStruct{1}{1}.length) =   pDif{idx}( pDif{idx}<-theoryStruct{1}{1}.length) +theoryStruct{1}{1}.length;

    [val{idx},sorti{idx}] = sort(posTrue);

end

figure,tiledlayout(2,1)
nm = {'(D)','(E)','(F)'}
for idx=1:length(pDif);
nexttile([1 1])
hold on
plot(val{idx},     pDif{idx}(sorti{idx}),'blackx')
    xlabel('Position (px)')
    ylabel('pd')

    title([nm{idx},'PD barcode island (',num2str(idx),')'], 'Interpreter','latex')

end
% Maybe use data from this for Figure 2 table?


% %
% f = figure
% tiledlayout(2,2,'TileSpacing','compact')
% nexttile
% hold on
% for i=1:length(barcodeGen)
%     plot(zscore(barcodeGen{i}.rawBarcode)+5*i,'black');
% end
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% title('(A) List of barcodes')
% 
% nexttile
% 
% plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1);
% title('(B) Overlap graph ')
% 
% nexttile([1 2])
% hold on
% 
% title('(C) Largest barcode island ', 'Interpreter','latex')
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
% 
% print('FIGS/Fig1.eps','-depsc','-r300');
% 
% % nexttile
% % hold on
% % 
% % title('(D) Barcode island 2 ')


%%
bb = 0.8;
    [newData] = update_sf_barset(barcodeIslandsData{idx}, bb,1);

        [newData2] = update_sf_barset(newData, 1/bb,1);

    figure,imagesc(newData2-barcodeIslandsData{idx});colorbar
    
end