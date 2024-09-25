function [] = final_figure_5(resRun,bG,theoryStruct)
%%
alphaNu1 = 0.09;
alphaNu2 = 0.085;
alphaN1 = 0.42;
alphaN2 = 0.004;
pthresh = 0.05;
nB = 175;
sets.minOverlap = 300;
sets.scDiffSetting = 0.05;
sets.pxDifSetting = 50;
sF = 0.8:0.025:1.2;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

ixtest = 8; % load the synthetic dataset with 3 islands




% maybe put into sep function
oS = resRun{ixtest}.oS;
barcodeGen = bG{ixtest}';
tS = theoryStruct{ixtest};

oS = oS(1:nB,1:nB);
barcodeGen = barcodeGen(1:nB);

[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu1, pthresh,[],alphaN1,alphaNu2,alphaN2); %localCCStruct

% % both on local and global
% [sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
%     calculate_sorted_pvals(oS,sets.minOverlap, alphaNu, pthresh); %localCCStruct

sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

% import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    Core.barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],0,3)
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

% import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   Core.create_graph_from_barsetgen(barsetGen,barcodeGen,[],0);


%% Run amplitude & iterative update
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
    [aaBarcodeIslandBlockRep{ix}, ~,~,clusterBarcodes{ix},~] = barcode_island_consensus(barcodeGen,cGenAll, ix, wminC,[],alphaNu1);

%     NN = 5;
    % wminC
    [matRep{ix},consensusNew{ix},barCon{ix},stats{ix},cIt{ix}] = iterative_readjustment_consensus(aaBarcodeIslandBlockRep{ix},barcodeGen, barIslands, cGenAll,ix,wminC,NN,alphaNu1);

end


%% Validation. 
sets.comparisonMethod = 'mass_pcc';
bT{1}.rawBitmask = ones(1,length(tS{1}.rawBarcode));
import Validation.scaled_pdif_vs_theory;
[val,sorti,valRes] = scaled_pdif_vs_theory(bT,tS{1},barcodeGen,cellfun(@(x) x{end}.idx,cIt,'UniformOutput',false),cellfun(@(x) x{end}.vals,cIt,'UniformOutput',false),sF,sets);
% import Validation.scaled_pdif_vs_theory;
% [val,sorti,valRes] = scaled_pdif_vs_theory(bT,tS{1},barcodeGen,barIslands,barcodeIslandsData,sF,sets);


% 
% 
% % tS.length = length(tS.rawBitmask);
% pDif = cell(1,length(barIslands));
% sorti =  cell(1,length(barIslands));
% val =  cell(1,length(barIslands));
% tablesb = cell(1,length(barIslands));
% for idx = 1:length(barIslands)
% 
%     %%1) Run comparison (all barcodes)
%     bars = barcodeGen(barIslands{idx}); 
%     barx = barIslands{idx};
%     sets.comparisonMethod = 'mass_pcc';
%     sets.nuF = 0.14;
%     [compI,rezI,~] = compare_to_t(bars,tS,sF,sets); % maybe calculate Stouffer score here also?
% 
%     [a] = cellfun(@(x) x.pval,compI);
%     [minV,minPos] = min(a);
%     %% 
%     selBarId = minPos;
% 
%     %% extract tables
%     bbS = compI{selBarId}.bestBarStretch;  % best bar stretch (bbS)
%     bbO = (compI{selBarId}.or(1)~=barcodeIslandsData{idx}(selBarId,3))+1;
% 
%     [bbSAll] = Core.update_sf_barset(barcodeIslandsData{idx}, bbS/barcodeIslandsData{idx}(selBarId,4), bbO);
% 
%     % need to also update positions for all the synCur barcodes based on
%     [tableS] = Core.synth_to_table(compI);
%     [tableSUpd] = Core.update_sf_barset(tableS, tableS(selBarId,4)/bbS, 1); % does nothing
% %     tableSUpd = tableS;
%     % update synth table. Want to have the same / instead of this, should
%     % update bbSall table?
%     % [bbSAll] = update_sf_barset(bbSAll, 1/tableS(selBarId,4), 1);
% 
% 
%     bbSAll(:,1:2) = bbSAll(:,1:2)-bbSAll(selBarId,1)+tableSUpd(selBarId,1); %
%    
%     bestPosFound = bbSAll(:,1);
%     
%     posTrue = tableSUpd(:,1);
% 
%     pDif{idx} =  min([sqrt((posTrue-bestPosFound).^2) sqrt((posTrue+theoryStruct{ixtest}{1}.length-bestPosFound).^2)...
%         sqrt((posTrue-theoryStruct{ixtest}{1}.length-bestPosFound).^2)]');
% 
%     [val{idx},sorti{idx}] = sort(posTrue);
%     tables{idx} = bbSAll;
% 
% end
% 
[cellfun(@(x) mean(x.pDif(x.pDif<30)),valRes)]
[cellfun(@(x) std(x.pDif(x.pDif<30)),valRes)]

cellfun(@(x) length(x.pDif(x.pDif>=30))/length(x.pDif),valRes)
cellfun(@(x) length(x.pDif(x.pDif>=300))/length(x.pDif),valRes)

%%
% idx = 1
% % figure
% % imagesc(min( tables{idx}(:,1)):max( tables{idx}(:,2)),1,theoryStruct.rawBarcode(min( tables{idx}(:,1)):max( tables{idx}(:,2))));colormap(gray)
% 
% import Core.consensus_from_table
% [outConsensus2] = consensus_from_table(bbSAll,bars);
% 
% 
% [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus2(x,:)),1,'first') +(find(~isnan(outConsensus2(x,:)),1,'last')-find(~isnan(outConsensus2(x,:)),1,'first'))/2,1:size(outConsensus2,1)));
% % 
% consensusToPlot1 = outConsensus2(idxv,:);
% 
% figure;imagesc(consensusToPlot1);colormap(gray)
% 
% 
% foundPos = min(bestPosFound);
% 
% shift = -min(bbSAll(:,1))+1;
% foundPosEnd = max(bbSAll(:,2));
% 
% subMat = consensusToPlot1(:,shift+foundPos:shift+foundPosEnd);
% subMat(all(isnan(subMat), 2), :) = [];
% 
% f = figure;
% tiledlayout(5,1)
% nexttile
% imagesc(tableS(selBarId,1):tableS(selBarId,2),1,theoryStruct{ixtest}{1}.rawBarcode(tableS(selBarId,1):tableS(selBarId,2)));colormap(gray)


%%

outConsensus = cellfun(@(x) x{end},matRep,'un',false);
f = figure;
% f=figure('Position', [10 10 714 500]);
tiledlayout(3,3,'TileSpacing','compact','Padding','none')



maxConsensus = max(cellfun(@(x) size(x,2),outConsensus));
% hold on

lettersAll = mat2cell('A':'Z',1,ones(1,26));

letters = lettersAll(1:end)
for idx=1:length(outConsensus)
    nexttile([1 3])
    hold on
    title(['(' letters{idx} ') Block representation (', num2str(idx),')'], 'Interpreter','latex')
        
    % sort based on starting position
    [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
    % 
    consensusToPlot1 = outConsensus{idx}(idxv,:);

    consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/2),1);
    consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
    imagesc([consensusBlock;consensusEmptyBlock;consensusToPlot1;]);colormap(gray)
    ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
%     text(size(consensusBlock,2)+10,size(consensusBlock,1)/2, 'Consensus','FontSize',10, 'FontName','Times',   'Color','black')
    axis off

    set(gca,'xtick',[])
    text(-650,size(consensusBlock,1)/2, 'Consensus','FontSize',7, 'FontName','Times',   'Color','black')

    xlim([1 maxConsensus])
end

nPixels = 1e5/500;% 500bp/px ;
x = [maxConsensus-nPixels-2020 maxConsensus-2020];
% x=[1 200]
y = [20.6 20.6 ]+30;
plot(x+100,y+30,'Linewidth',2,'Color','white')
% text(x(1),y(1)+0.1, '10 $\mu m$','FontSize',12, 'Color','white','Interpreter','latex')
text(x(1)+nPixels/10,y(1)+40, '100 kb','FontSize',8, 'FontName','Times',   'Color','white')
set(gcf, 'Color', 'w')




%
% 
% nexttile([1 3])
% 
% maxConsensus = max(cellfun(@(x) size(x,2),outConsensus));
% hold on
% 
% title('(A) Block representation for barcode island (1)', 'Interpreter','latex')
% 
% idx = 1;
% 
% % sort based on starting position
% [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% % 
% consensusToPlot1 = outConsensus{idx}(idxv,:);
% 
% % imagesc(consensusToPlot1);colormap(gray)
% % xlim([1 size(consensusToPlot1,2)])
% % ylim([1 size(consensusToPlot1,1)])
% 
%     consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/5),1);
%     consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
% %     imagesc([consensusToPlot1;consensusEmptyBlock;consensusBlock]);colormap(gray)
%     imagesc([consensusBlock;consensusEmptyBlock;consensusToPlot1;]);colormap(gray)
% 
% %     set(gca,'YDir','normal')
%     ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
% 
% 
% % imagesc(outConsensus{idx}(idxv,:));colormap(gray)
% axis off
% 
% set(gca,'xtick',[])
% text(-500,size(consensusBlock,1)/2, 'Consensus','FontSize',7, 'FontName','Times',   'Color','black')
% 
% xlim([1 maxConsensus])
% 
% nexttile([1 3])
% hold on
% idx = 2;
% 
% % sort based on starting position
% [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% % 
% consensusToPlot2 = outConsensus{idx}(idxv,:);
% 
% % imagesc(consensusToPlot2);colormap(gray)
% xlim([1 maxConsensus])
% ylim([1 size(consensusToPlot2,1)])
% % 
% % imagesc(outConsensus{idx}(idxv,:));colormap(gray)
% consensusToPlot1 = outConsensus{idx}(idxv,:);
% 
%   consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/5),1);
%     consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
% %     imagesc([consensusToPlot1;consensusEmptyBlock;consensusBlock]);colormap(gray)
%     imagesc([consensusBlock;consensusEmptyBlock;consensusToPlot1;]);colormap(gray)
% 
% %     set(gca,'YDir','normal')
%     ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
% 
% 
% 
% axis off
% 
% % set(gca,'xtick',[])
% 
% % xlim([1 size(consensusToPlot2,2)])
% 
% title('(B) Block representation for barcode island (2)', 'Interpreter','latex')
% 
% text(-500,size(consensusBlock,1)/2, 'Consensus','FontSize',7, 'FontName','Times',   'Color','black')
% 
% nexttile([1 3])
% hold on
% idx = 3;
% 
% % sort based on starting position
% [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
% % 
% consensusToPlot3 = outConsensus{idx}(idxv,:);
% 
% % imagesc(consensusToPlot3);colormap(gray)
% % xlim([1 size(consensusToPlot3,2)])
% ylim([1 size(consensusToPlot3,1)])
% xlim([1 maxConsensus])
% 
% % imagesc(outConsensus{idx}(idxv,:));colormap(gray)
% consensusToPlot1 = outConsensus{idx}(idxv,:);
% 
%   consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/5),1);
%     consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
% %     imagesc([consensusToPlot1;consensusEmptyBlock;consensusBlock]);colormap(gray)
%     imagesc([consensusBlock;consensusEmptyBlock;consensusToPlot1;]);colormap(gray)
% 
% %     set(gca,'YDir','normal')
%     ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
% 
% 
% 
% % set(gca,'xtick',[])
% set(gca,'ytick',[])
% axis off
% 
% % set(gca,'ytick',[])
% xlim([1 maxConsensus])
% 
% xlabel('Position (kb)')
% 
% xticksOriginal = get(gca, 'XTick');
% xtickLabelsOriginal = get(gca, 'XTickLabel');% nexttile([1 2])
% 
% scaleFactor = 0.5;
% % Scale the tick locations
% xticksScaled = xticksOriginal * scaleFactor;
% 
% xtickLabelsScaled = cellstr(num2str(xticksScaled'));
% set(gca, 'XTickLabel', xtickLabelsScaled);
% % 
% text(-500,size(consensusBlock,1)/2, 'Consensus','FontSize',7, 'FontName','Times',   'Color','black')
% 
% % xlim([1 size(consensusToPlot3,2)])
% 
% title('(C) Block representation for barcode island (3)', 'Interpreter','latex')
% 
% 
% 
% % nm = {'(D)','(E)','(F)'}
% % for idx=1:3

    nexttile([1 3])

    par = {'x','o','|'      };% nexttile
hold on
for idx=1:length(val)
    

    plot(valRes{idx}.posTrue,valRes{idx}.bbSAll(:,1),par{idx});xlabel('Ground truth position (px)','Interpreter','latex');ylabel('Placement (px)','Interpreter','latex');

%     plot(val{idx},     pDif{idx}(sorti{idx}),'blackx')
%     xlabel('Position (px)')
%     ylabel('pd')

%     title([nm{idx} ' Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')
    xlim([1 length(theoryStruct{ixtest}{1}.rawBarcode)])
    title('(C) Dot plot for barcode islands','Interpreter','latex');
end
legend({'(1)','(2)','(3)'},'Interpreter','latex','Location','best')

%%

% print('-depsc','-tiff','-r350', '-vector','FIGS/Fig5.eps')%large filesize
% exportgraphics(gcf,,'BackgroundColor','none','ContentType','vector')
    print('FIGS/Fig5.eps','-depsc','-r600');


%% Extra
% 
idx = 2;
figure
hold on
plot(val{idx},     valRes{idx}.pDif(sorti{idx}),'blackx')
    xlabel('Position (px)')
    ylabel('pd')

    title(['Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')

%%

% figure
% histogram(valRes{1}.pDif);
% % hold on
% histogram(valRes{2}.pDif);
% % histogram(valRes{3}.pDif);


%%
% print('FIGS/Fig5.eps','-depsc','-r300');

%% Suplementary figure S6

import Plot.ref_based_assembly_islands;
[f] = ref_based_assembly_islands(barcodeGen,barIslands,synthStr{ixtest},theoryStruct{ixtest})

print(['FIGS/FigS6','.eps'],'-depsc','-r500');

%% 
for idx = 1:2;
  bars = barcodeGen(barIslands{idx}); 
    barx = barIslands{idx};
   
[outConsensus2, coverage2, pval3,f] = gen_reference_based_assembly(bars,synthStr{ixtest}(barx),theoryStruct{ixtest},'test11',inf);
    print(['FIGS/Fig5S_', num2str(idx), '.eps'],'-depsc','-r500');

end




%% What happens at the edge of two islands?

positions = [3000 3000]; % positions between the two that we are interested in, do they match up or not ?

featureLen  = 3000;
f=figure('Position', [10 10 714 500]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')
import Validation.plot_compare_to_thry;
 plot_compare_to_thry(f,positions(1),1,valRes,cIt,barcodeGen,barcodeIslandsData,wminC,tS,featureLen,1)

 import Validation.plot_compare_to_thry;
 plot_compare_to_thry(f,positions(2),2,valRes,cIt,barcodeGen,barcodeIslandsData,wminC,tS,featureLen,0)
     print('FIGS/FigS20.eps','-depsc','-r500');

