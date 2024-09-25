function [] = plot_compare_to_thry(f,foundPos,ix,valRes,cIt,barcodeGen,barcodeIslandsData,wminC,tS,featureLen,plotT)

% find all barcodes in compI matching the required position (foundPos)
allPos = cellfun(@(x) x.pos(1), valRes{ix}.compI);
allPosCoef = cellfun(@(x) x.maxcoef(1), valRes{ix}.compI);

closeregion = find(and(allPos>foundPos,foundPos<foundPos+300));

[score,selId] = max(allPosCoef(closeregion));

selBar = closeregion(selId);

import Validation.scaled_specific_barcode_vs_theory;
[valResSpec{ix},val,sorti] = scaled_specific_barcode_vs_theory(valRes{ix}.compI,selBar,barcodeIslandsData,tS{1},ix);

selBarId = valResSpec{ix}.selBarId;
bbSAll = valResSpec{ix}.bbSAll;
tableS = valResSpec{ix}.tableS;
posTrue =  valResSpec{ix}.posTrue;

% create barcode island consensuses for all barcode islands from the best
% positions / single iteration (since should already be corrected, we just
% slightly re-scaled according to theory)
import Plot.islandsPosStruct;
import Core.barcode_island_consensus;
pS = cell(1,length(valResSpec));
% for i=1:length(valResSpec)
    pS{ix}.comparisonStruct = islandsPosStruct({valResSpec{ix}.bbSAll},{cIt{ix}{end}.idx});
    pS{ix}.idx = cIt{ix}{end}.idx;
    [aaRep{ix}, ~ , ~,clusterBarcodesaa{ix},rawBarcodeIslandBlockRep] = barcode_island_consensus(barcodeGen,pS, ix, wminC);
% end
% 
outConsensus2 = aaRep{ix};

[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus2(x,:)),1,'first') +(find(~isnan(outConsensus2(x,:)),1,'last')-find(~isnan(outConsensus2(x,:)),1,'first'))/2,1:size(outConsensus2,1)));
% 
consensusToPlot2 = outConsensus2(idxv,:);
%

if bbSAll(selBar,1)+foundPos> tS{1}.length
    bbSAll(:,1:2) =  bbSAll(:,1:2)-tS{1}.length;
end
startSelBar = bbSAll(selBar,1)-min(bbSAll(:,1))+1;
starZone =  startSelBar + foundPos-bbSAll(selBar,1);
stopZoneBar = min(starZone+featureLen-1,size(consensusToPlot2,2));
starZoneBar = max(1,starZone);

% stopZone = stopBar+ foundPos-bbSAll(selBar,1);
% featStop = featStart;

% shift = -min(bbSAll(:,1))+1;

% indexes = [max(1,starZone) min(size(consensusToPlot2,2),shift+foundPos+bbSAll(selBarId,2)-bbSAll(selBarId,1))];
indexes = [starZoneBar stopZoneBar];

if starZoneBar~=starZone
    start = -starZone+1;
else
    start = 1;
end

if  stopZoneBar ~= starZone+featureLen-1
    stop = indexes(2)-indexes(1)+1;
else 
    stop = start+stopZoneBar-1;
end
 indexes2 = foundPos+[start stop]-1;

% indexes = indexes+shift000;
subMat = consensusToPlot2(:,indexes(1):indexes(2));
subMat(all(isnan(subMat), 2), :) = [];
%
bar = [tS{1}.rawBarcode tS{1}.rawBarcode];
% tableS= bbSAll;
if plotT
nexttile
% imagesc(indexes2(1):indexes2(2),1,bar(indexes2(1):indexes2(2)));colormap(gray)
imagesc(foundPos:foundPos+featureLen-1,1,bar(foundPos:foundPos+featureLen-1));colormap(gray)

% imagesc(tableS(selBarId,1):tableS(selBarId,2),1,bar(tableS(selBarId,1):tableS(selBarId,2)));colormap(gray)
axis off
title(['(A) Part of theory corresponding to ', num2str(foundPos),'-',num2str(foundPos+featureLen-1), ' (px)'],'Interpreter','latex')
% xlim([1 12000])
% xlim([indexes2(1)+start indexes2(2)])
%
xlim([foundPos foundPos+featureLen-1])
end

ax1=nexttile;
% subMat(subMat<-4) = nan; % mask edge cases for better visibility
% subMat(subMat>4) = nan;
hold on
xlim([foundPos foundPos+featureLen-1])

% zoomin
consensusToPlot1 = subMat;
consensusToPlot1 = consensusToPlot1(find(sum(~isnan(consensusToPlot1'))),:);


consensusBlock  = repmat(mean(consensusToPlot1,'omitnan'),round(size(consensusToPlot1,1)/5),1);
consensusEmptyBlock  = max(consensusToPlot1(:)).*ones(round(size(consensusBlock,1)/5),size(consensusBlock,2));
%     imagesc([consensusToPlot1;consensusEmptyBlock;consensusBlock]);colormap(gray)
imagesc(indexes2(1):indexes2(2),1,[consensusBlock; consensusEmptyBlock;consensusToPlot1]);colormap(gray)


% imagesc(tableS(selBarId,1):tableS(selBarId,2),1,subMat);colormap gray

    ylim([0.5 size([consensusBlock;consensusEmptyBlock;consensusToPlot1;],1)+0.5])
%     text(tableS(selBarId,1)-100,size(consensusBlock,1)/2, 'Consensus','FontSize',9, 'FontName','Times',   'Color','black')
%     xlim([foundPos:foundPos+1000-1])

title(['Corresponding part of block representation of barcode island (',num2str(ix),')'],'Interpreter','latex')
xlabel('Pos (px) along theory','Interpreter','latex')
ax1.YAxis.Visible = 'off'; 
axis off;
% nexttile
% xlim([indexes2(1) indexes2(2)])


end

