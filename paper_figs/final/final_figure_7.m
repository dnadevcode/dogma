% validation of fig_6


% First generate theory:
w = 300;% minimum overlap length
sF = 0.8:0.025:1.2;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

fastas = {'018_final_polish.fasta','DA32087.fasta','DA68335.fasta'};
theoryIdxs = {1, 3, 2, nan, 2 , 2, 3}; % known theory indexes
dataSetIdx = 3; % 5,... test datasets (part of the data)

nmPerPx = 110;
nmbp = 0.25; % ? 

% quickly calculate theory for given single fasta file
[theoryStructRev,theoryStruct,bT] = prep_thry_for_local_comp(fastas(theoryIdxs{dataSetIdx}), nmbp, nmPerPx, 1);

% Compare to theory
sets.comparisonMethod = 'mass_pcc';
sets.nuF = 0.085;
sets.nF = 0.15;
%%
% Based on best scoring position. Can change to based on
% consensus/consensus average etc?
import Validation.scaled_pdif_vs_theory;
% [val,sorti,valRes] = scaled_pdif_vs_theory(bT,theoryStruct,barcodeGen,barIslands,barcodeIslandsData,sF,sets);
[val,sorti,valRes] = scaled_pdif_vs_theory(bT,theoryStruct,barcodeGen,cellfun(@(x) x{end}.idx,cIt,'UniformOutput',false),cellfun(@(x) x{end}.vals,cIt,'UniformOutput',false),sF,sets);


% 
[cellfun(@(x) mean(x.pDif(x.pDif<30)),valRes)]
[cellfun(@(x) std(x.pDif(x.pDif<30)),valRes)]

cellfun(@(x) length(x.pDif(x.pDif>=500))/length(x.pDif),valRes)

%% non-rescaled
% [valNonA,sortiNonA,valResNonA] = scaled_pdif_vs_theory(bT,theoryStruct,barcodeGen,barIslands,barcodeIslandsData,sF,sets);
% % 
% [cellfun(@(x) mean(x.pDif(x.pDif<30)),valResNonA)]
% [cellfun(@(x) std(x.pDif(x.pDif<30)),valResNonA)]
% 
% cellfun(@(x) length(x.pDif(x.pDif>=300))/length(x.pDif),valResNonA)


%%
ix = 1;
selBarId = valRes{ix}.selBarId;
bbSAll = valRes{ix}.bbSAll;
tableS = valRes{ix}.tableS;
posTrue =  valRes{ix}.posTrue;


% create barcode island consensuses for all barcode islands from the best
% positions / single iteration (since should already be corrected, we just
% slightly re-scaled according to theory)
import Plot.islandsPosStruct;
import Core.barcode_island_consensus;
pS = cell(1,length(valRes));
for i=1:length(valRes)
    pS{i}.comparisonStruct = islandsPosStruct({valRes{i}.bbSAll},{cIt{i}{end}.idx});
    pS{i}.idx = cIt{i}{end}.idx;
    [aaRep{i}, ~ , ~,clusterBarcodesaa{i},rawBarcodeIslandBlockRep] = barcode_island_consensus(barcodeGen,pS, i, wminC);

end

% [aaBarcodeIslandBlockRep, Z , allscores,clusterBarcodes,rawBarcodeIslandBlockRep] = barcode_island_consensus(barcodeGen,pS, ix, wminC);
% outConsensus2 = [];
% for i=1:length(valRes)
outConsensus2 = aaRep{ix};
% end
% scaled table to map:
% import Core.consensus_from_table; % todo: better scaling
% [outConsensus2] = consensus_from_table(valRes{ix}.bbSAll,valRes{ix}.bars);
% [outConsensus2] = consensus_from_table(barcodeIslandsData{2},valRes{ix}.bars);

[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus2(x,:)),1,'first') +(find(~isnan(outConsensus2(x,:)),1,'last')-find(~isnan(outConsensus2(x,:)),1,'first'))/2,1:size(outConsensus2,1)));
% 
consensusToPlot2 = outConsensus2(idxv,:);
% consensusToPlot2(consensusToPlot2<0) = nan;
% figure
% imagesc(consensusToPlot2);colormap(gray)
% figure;imagesc(consensusToPlot1);colormap(gray)


foundPos = valRes{ix}.bestPosFound(selBarId);

shift = -min(valRes{ix}.bbSAll(:,1))+1;

subMat = consensusToPlot2(:,shift+foundPos:shift+foundPos+bbSAll(selBarId,2)-bbSAll(selBarId,1));
subMat(all(isnan(subMat), 2), :) = [];
%%
tableS= bbSAll;
f=figure('Position', [10 10 714 500]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')
nexttile
imagesc(tableS(selBarId,1):tableS(selBarId,2),1,theoryStruct.rawBarcode(tableS(selBarId,1):tableS(selBarId,2)));colormap(gray)
axis off
title('(A) Part of theory corresponding to optimally matching barcode of barcode island (1)','Interpreter','latex')
% xlim([1 12000])

ax1=nexttile;
% subMat(subMat<-4) = nan; % mask edge cases for better visibility
% subMat(subMat>4) = nan;
hold on

    % zoomin
    consensusToPlot1 = subMat;
    consensusToPlot1 = consensusToPlot1(find(sum(~isnan(consensusToPlot1'))),:);
    

    consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/5),1);
    consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
%     imagesc([consensusToPlot1;consensusEmptyBlock;consensusBlock]);colormap(gray)
    imagesc(tableS(selBarId,1):tableS(selBarId,2),1,[consensusBlock; consensusEmptyBlock;consensusToPlot1]);colormap(gray)


% imagesc(tableS(selBarId,1):tableS(selBarId,2),1,subMat);colormap gray

    ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
    text(tableS(selBarId,1)-110,size(consensusBlock,1)/2, 'Consensus','FontSize',7, 'FontName','Times',   'Color','black')
xlim([tableS(selBarId,1) tableS(selBarId,2)])

title('(B) Corresponding part of block representation of barcode island (1)','Interpreter','latex')
xlabel('Pos (px) along theory','Interpreter','latex')
ax1.YAxis.Visible = 'off'; 
nm = {'(C)','(D)','(E)','(F)'};

% xlim([1 12000])
nexttile
hold on
for idx=1:length(val)
    

    plot(valRes{idx}.posTrue,valRes{idx}.bbSAll(:,1),'x');xlabel('Ground truth position (px)','Interpreter','latex');ylabel('Placement (px)');

%     plot(val{idx},     pDif{idx}(sorti{idx}),'blackx')
%     xlabel('Position (px)')
%     ylabel('pd')

%     title([nm{idx} ' Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')
    xlim([1 length(theoryStruct.rawBarcode)])
end
    title('(C) Dot plot for barcode islands','Interpreter','latex');

legend({'(1)','(2)','(3)'},'Interpreter','latex')
% if bbSAll
    print('FIGS/Fig7.eps','-depsc','-r300');

    %% Specific part
    positions = [2000]; % positions between the two that we are interested in, do they match up or not ?

featureLen  = 3500;
f=figure('Position', [10 10 714 500]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact')
import Validation.plot_compare_to_thry;
ii =1;
 plot_compare_to_thry(f,positions(1),ii,valRes,cIt,barcodeGen,cellfun(@(x) x{end}.vals,cIt,'UniformOutput',false),wminC,{theoryStruct},featureLen,1)
ii =2;
 import Validation.plot_compare_to_thry;
 plot_compare_to_thry(f,positions(1),ii,valRes,cIt,barcodeGen,cellfun(@(x) x{end}.vals,cIt,'UniformOutput',false),wminC,{theoryStruct},featureLen,0)
     print('FIGS/Fig7S_1.eps','-depsc','-r500');
%%
  
    %% S10
%     idx = 2;

f=figure('Position', [10 10 714 200]);

plot(tableS(selBarId,1):tableS(selBarId,2), zscore(nanmean(subMat)));hold on;plot(tableS(selBarId,1):tableS(selBarId,2),zscore(theoryStruct.rawBarcode(tableS(selBarId,1):tableS(selBarId,2))))
xlabel('Pos (px) along theory','Interpreter','latex')
legend({'Consensus for barcode island 1','Theory for barcode island 1'},'Interpreter','latex')
print('FIGS/FigS10.eps','-depsc','-r500');


% f=figure('Position', [10 10 714 200]);
% tiledlayout(size(aaRep,1),1);
% for i=1:size(aaRep,1)
%     nexttile
%     outConsensus2 = aaRep{ix};
%     [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus2(x,:)),1,'first') +(find(~isnan(outConsensus2(x,:)),1,'last')-find(~isnan(outConsensus2(x,:)),1,'first'))/2,1:size(outConsensus2,1)));
%     % 
%     consensusToPlot2 =  (nanmean(outConsensus2(idxv,:)));
% %     plot(tableS(selBarId,1):tableS(selBarId,2), zscore(nanmean(subMat)));hold on;plot(tableS(selBarId,1):tableS(selBarId,2),zscore(theoryStruct.rawBarcode(tableS(selBarId,1):tableS(selBarId,2))))
%     
%     selBarId = valRes{i}.selBarId;
%     bbSAll = valRes{i}.bbSAll;
%     tableS = valRes{i}.tableS;
% 
% %     plotdiff = bbSAll(selBarId,1)
% 
%     plot((consensusToPlot2-nanmean(consensusToPlot2))./nanstd(consensusToPlot2))
%     hold on
%     plot(zscore(theoryStruct.rawBarcode))
% %     plot(tableS(selBarId,1):tableS(selBarId,2), zscore(nanmean(subMat)));hold on;plot(tableS(selBarId,1):tableS(selBarId,2),zscore(theoryStruct.rawBarcode(tableS(selBarId,1):tableS(selBarId,2))))
% 
% 
% end
%%
ix = 1;
selBarId = valRes{ix}.selBarId;
bbSAll = valRes{ix}.bbSAll;
tableS = valRes{ix}.tableS;

    %% S8
    ix = 2;
import Core.barcode_island_consensus;
[aaBarcodeIslandBlockRep, Z , allscores,clusterBarcodes,rawBarcodeIslandBlockRep] = barcode_island_consensus(barcodeGen,cGenAll, ix, wminC);

% run some extra steps
NN = 10;
% wminC
[matRep,consensusNew,barCon,stats,cGen] = iterative_readjustment_consensus(aaBarcodeIslandBlockRep,barcodeGen, barIslands, cGenAll,ix,wminC,NN)

f = consensus_ci(aaBarcodeIslandBlockRep,rawBarcodeIslandBlockRep,matRep{end})
print('FIGS/FigS9.eps','-depsc','-r400');

% consensus_ci(aaBarcodeIslandBlockRep,rawBarcodeIslandBlockRep)
%%
idx = 2;
goodBars = logical(cellfun(@(x) x.pval,valRes{idx}.compI,'un',true)<0.001);
badBars = ~goodBars;
% figure,plot(valRes{idx}.pDif(sorti{idx}(goodBars)))
% figure,plot(valRes{idx}.pDif(sorti{idx}(badBars)))

%%
figure,plot(posTrue,bbSAll(:,1),'x');xlabel('ground truth (px)');ylabel('Placement (px)');title('Dot plot for barcode island (2)')

% Local similarity
wSim = 500;
[MPthr, MPIthr1] = local_compare_thry(wSim,[], nmbp, nmPerPx, 1,32,theoryStructRev)
  
mpVals = MPthr{1}(1:theoryStruct.length-wSim+1);
% mpVals = MPthr{1}(1:theoryStruct.length);


figure
tiledlayout(2,1)
nexttile
idx =2
plot(val{idx},     pDif{idx}(sorti{idx}),'blackx')
    xlabel('Position (px)')
    ylabel('pd')

    title([' Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')
    xlim([1 length(theoryStruct.rawBarcode)])

    nexttile
    hold on
title('Local similarity along the theory barcode')
plot(mpVals)
    xlabel('Position (px)')

    xlim([1 length(theoryStruct.rawBarcode)])


    %%
    % 
idx = 1;
figure
hold on
plot(val{idx},     valRes{idx}.pDif(sorti{idx}),'blackx')
    xlabel('Position (px)')
    ylabel('pd')

    title(['Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')