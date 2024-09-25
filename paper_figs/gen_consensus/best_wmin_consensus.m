% Test to see what's the best wmin for getting barcode_island_consensus


ix = 1;
wvals = 10:10:200;
wmin = 300;
sets.comparisonMethod = 'mpnan';
sets.nuF = 0.08;
sets.w = 300;
NN = length(wvals);

sF = 0.85:0.01:1.15;
pccMean = zeros(1,NN);

for i = 1:NN
    wminC = wvals(i);

import Core.barcode_island_consensus;
[scoreAlignedData, Z , allscores] = barcode_island_consensus(barcodeGen,cGenAll, ix, wminC);

% Temporaro consensus using nanmean
con1 = mean(scoreAlignedData,'omitnan');
% conz = mean(outConsensus{ix},'omitnan');
% figure,plot((con1-nanmean(con1))/nanstd(con1,1))
% hold on,plot((conz-nanmean(conz))/nanstd(conz,1))
% legend({'New scaling','Z-score'}) 
% Compare the consensus barcode against theory using matrix profile (MP) 

bar =[];
bar{1}.rawBarcode = con1;
bar{1}.rawBitmask = ~isnan(con1);
% bar{2}.rawBarcode = conz;
% bar{2}.rawBitmask = ~isnan(conz);

% this only gets the best coefficient
% [compI,rezI,~] = compare_to_t(bar,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?
% 
% import Validation.scaled_pdif_vs_theory;
% [val,sorti,valRes] = scaled_pdif_vs_theory(bT,theoryStruct,barcodeGen,barIslands,barcodeIslandsData,sF,sets);
import Core.compare_rescaled_theory;
% compare_rescaled_theory(theoryStruct,sF, wmin,numWorkers)
[mp1, mpI, indexesT, indexesT2] = compare_rescaled_theory(bar, theoryStruct,sF, wmin,numWorkers);

% figure,
% tiledlayout(2,1)
% nexttile
% plot(mp1(indexesT{1}(1):indexesT{1}(2)));hold on;
% plot(mp1(indexesT{2}(1):indexesT{2}(2)));
% ylabel('PCC')
% legend({'New consensus','Old consensus'},'Location','eastoutside')
% nexttile
% plot(sum(~isnan(scoreAlignedData)));
% legend({'Coverage'},'Location','eastoutside')
    pccMean(i)  = mean(mp1(indexesT{1}(1):indexesT{1}(2)),'omitnan');
pccMean(i)
end

figure,plot(wvals,pccMean)
xlabel('w')
ylabel('Mean pcc');

