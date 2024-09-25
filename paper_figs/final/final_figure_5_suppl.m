function [] = final_figure_5_suppl(resRun,bG,theoryStruct, nB,ixtest,imnr)

if nargin < 4
    nB = 175;
    ixtest = 7;
    imnr = 21;
end
%%
alphaNu1 = 0.09;
alphaNu2 = 0.085;
alphaN1 = 0.42;
alphaN2 = 0.004;
pthresh = 0.05;
% nB = 175; % the only diff
sets.minOverlap = 300;
sets.scDiffSetting = 0.05;
sets.pxDifSetting = 50;
sF = 0.9:0.025:1.1;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

% ixtest = 8; % load the synthetic dataset with 3 islands




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

% 
[cellfun(@(x) mean(x.pDif(x.pDif<20)),valRes)]
[cellfun(@(x) std(x.pDif(x.pDif<20)),valRes)]

cellfun(@(x) length(x.pDif(x.pDif>=20))/length(x.pDif),valRes)
cellfun(@(x) length(x.pDif(x.pDif>=300))/length(x.pDif),valRes)

%%

outConsensus = cellfun(@(x) x{end},matRep,'un',false);
f = figure;
tiledlayout(1+size(outConsensus,2),3,'TileSpacing','compact','Padding','none')



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
    title(['(', letters{length(outConsensus)+1} , ') ',' Dot plot for barcode islands'],'Interpreter','latex');
end
legend({'(1)','(2)','(3)'},'Interpreter','latex','Location','best')

%%
print(['FIGS/Fig' ,num2str(imnr),'.eps'],'-depsc','-r500');

%% Extra
% 
% idx = 1;
% figure
% hold on
% plot(val{idx},     valRes{idx}.pDif(sorti{idx}),'blackx')
%     xlabel('Position (px)')
%     ylabel('pd')
% 
%     title(['Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')

%%
%%
% print('FIGS/Fig5.eps','-depsc','-r300');

%% Suplementary figure S6

import Plot.ref_based_assembly_islands;
[f] = ref_based_assembly_islands(barcodeGen,barIslands,synthStr{ixtest},theoryStruct{ixtest})

print(['FIGS/FigS6','.eps'],'-depsc','-r500');

%% 
% for idx = 1:2;
%   bars = barcodeGen(barIslands{idx}); 
%     barx = barIslands{idx};
%    
% [outConsensus2, coverage2, pval3,f] = gen_reference_based_assembly(bars,synthStr{ixtest}(barx),theoryStruct{ixtest},'test11',inf);
%     print(['FIGS/Fig5S_', num2str(idx), '.eps'],'-depsc','-r500');
% 
% end




%% What happens at the edge of two islands?

positions = [8000]; % positions between the two that we are interested in, do they match up or not ?

featureLen  = 2000;
f=figure('Position', [10 10 714 500]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')
import Validation.plot_compare_to_thry;
 plot_compare_to_thry(f,positions(1),1,valRes,cIt,barcodeGen,barcodeIslandsData,wminC,tS,featureLen,1)

 import Validation.plot_compare_to_thry;
 plot_compare_to_thry(f,positions(2),2,valRes,cIt,barcodeGen,barcodeIslandsData,wminC,tS,featureLen,0)
     print('FIGS/FigS22.eps','-depsc','-r500');



