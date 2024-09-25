% evaluation figure.
%
% Individual barcodes can be compare against the theory, using
%
% 1) Standard PCC (+ pvalue to sort)
% 2) MP (+ p-value to sort)
% 3) Full discriminative barcode analysis against whole database
%
% As output, gives a position and significance


% For fig_evaluation, we compare the theory barcode against all the
% experimental barcodes at once (don't separate them into different
% folders). This allows us to use a single nm/bp length re-scaling factor.


% psffac = 1; % scaling factor (in case PSF needs to be something else than 300nm)
numWorkers = 30; % num workers, in one node there are 30
minLen = [150:50:3000]; % min to max % could take more points for more accurate..
sets.comparisonMethod = 'mass_pcc';
sF = 0.8:0.01:1.2;
% nmPerPx = 110;

% choose null-model
nullModel='pval'; % alternative 'mpmax', but we stick to pval as default

% nmbp = mean(cellfun(@(x)   x{1}.nmBpidFold,kymoStructs));

nmbp = 0.25;
import Thry.gen_theoretical;
[thrInd,~,~] = gen_theoretical(fastas(theoryIdxs{dataSetIdx}),nmbp,0,kymoStructs{1}{1}.nmpxidFold); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly

import Zeromodel.calc_null_model_thresh;
  [mpMaxLenBasedInd] = calc_null_model_thresh(nmbp,fastas,...
    theoryIdxs,dataSetIdx,nullModel,kymoStructs{1}{1}.nmpxidFold,numWorkers,minLen,{thrInd},1,0);

% compare to theory: either pcc or mp
[comparisonStructInd,~,~] = compare_to_t(barcodeGen,thrInd,sF,sets);

% comparison struct to overlap struct:
import Core.comparisonsres_to_struct; 
% [oSreal] = comparisonsres_to_struct(barcodeGen, comparisonStructInd,thrInd.length);


allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStructInd);
allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStructInd);

idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
coefDiffs = allCoefs-mpMaxLenBasedInd{1}(idxThresh);

passthreshInd = find(coefDiffs > 0);
successRateInd(1) = length(passthreshInd)/length(coefDiffs); % this could still have false positives.. ?

% alternative
pasthresh = find(cellfun(@(x) x.pval,comparisonStructInd) < 0.0001);

if plotExample
    import CBT.Hca.UI.Helper.plot_any_bar;
    plot_any_bar(4,barcodeGen,{comparisonStructInd},thrInd,1);
end



bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStructInd(passthreshInd), [], [], [],cumsum(cellfun(@(x) x.length,{thrInd})),1);%1/bpPx*10^6

bg = barcodeGen(passthreshInd);
[outConsensus, coverage, pval] = gen_reference_based_assembly(bg,comparisonStructInd(passthreshInd),{thrInd},'test11',inf);

[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,comparisonStructInd,{thrInd},'test11',inf);


%%

% Generate CB theory 
% nmpx = 110;
nmPerPx = 110; % nm/px (can be extracted from image info file)
nmbp = 0.225; % nm/bp (can be extracted via dna_barcode_matchmaker lambda pipeline
psffac = 1; % scaling factor (in case PSF needs to be something else than 300nm)
numWorkers = 30; % num workers, in one node there are 30
lenSeq = 0;
fastaFile =  {'DA32087.fasta'};
fastaFile = {'018_final_polish.fasta'};
% fastaFile = {'GCF_002076835.1_ASM207683v1_genomic.fna'};
addpath('/export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/')

% [theoryStructRev,theoryStruct,barcodeGen] = prep_thry_for_local_comp(fastaFile, nmbp, nmPerPx, psffac);



%% Minimum length plot (possibly as heatmap?) 
minLen =[200:50:3000]; % min to max
[MP,mpMaxLenBased,theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,psffac,numWorkers,minLen, lenSeq);

data = importdata(theoryStructRev{1}.filename);

ii = 1;
figure,findpeaks(MP{ii}{1}(1:end/2),'MinPeakDistance',200,'SortStr','descend')
xlabel('Position (genome)')
ylabel('PCC score')

[a,b] = findpeaks(MP{ii}{1},'MinPeakDistance',200,'SortStr','descend');

% [a,b] = sort(MP{1}{1},'desc','MissingPlacement','last')
f = figure,plot(minLen,mpMaxLenBased)
xlabel('Overlap length','Interpreter','latex');
ylabel('$C_{thresh}$','Interpreter','latex')

%%
nmbpvals = 0.15:0.01:0.3;
mpMax = zeros(length(nmbpvals),length(minLen));
for ii=1:length(nmbpvals)
    [MP,mpMax(ii,:),theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbpvals(ii),psffac,numWorkers,minLen, lenSeq);
end
f = figure;
heatmap(minLen,nmbpvals,mpMax)
title('Heatmap for the local similarities of barcodes')
% xlabel('Overlap length','Interpreter','latex');
% ylabel('nm/bp extension factor','Interpreter','latex')

print('output/FigS1.eps','-depsc','-r300');

%% Individual map: possibly as supplementary figure
ix = 1;
pos2 = MPI{ii}{1}(b(ix)+1);

bar1 = data(b(ix):b(ix)+minLen(ii)-1);
bar2 =  data(pos2:pos2+minLen(ii)-1);

lastPx = find(isnan(data));

figure,plot(bar1)
hold on
plot(bar2)
legend({['Starts ' num2str(b(ix))],['Starts ' num2str(pos2)]},'location','southoutside')
title(['PCC = ' num2str(a(ix) )])

pcc = @(x,y) zscore(x,1)*zscore(y,1)'/length(x)

pcc(bar1,bar2)

%% Supplemnetary: BLAST:
% /export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/

%% Local similarity



% nmbp = 0.24;
% psffac = 1;
% nmbp = nmbp*psffac;
% nmpx = nmpx*psffac;

import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmPerPx);

sets.comparisonMethod = 'mass_pcc';

sF = 0.9:0.01:1.1;
[comparisonStruct,rezMax,bestBarStretch] = compare_to_t(bgAll,theoryStruct,sF,sets);

allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);

idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
coefDiffs = allCoefs-mpMaxLenBased(idxThresh);

barsPassThresh = find(coefDiffs > 0.05);
% barsPassThresh = find(ones(1,length(comparisonStruct)));

super_quick_plot(6,bgAll,comparisonStruct,theoryStruct)

% 
% import Zeromodel.beta_ev_pdf;
%  p = beta_ev_pdf(xx, a_fit2(k), 1, n_fit2(k));%y_est(k)

bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct(barsPassThresh), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

bg = bgAll(barsPassThresh);
[outConsensus, coverage, pval] = gen_reference_based_assembly(bg,comparisonStruct(barsPassThresh),theoryStruct,'test11',inf);


% quick_visual_plot(46,1,bgAll,rezMax,bestBarStretch,theoryStruct)

[compStr,~,calcLengths] = compare_to_t_mp(bgAll,theoryStruct,sF,300); %todo: extend local to full
[comparisonStruct{4}.pos(1) compStr{2}.pos(1)]

passthreshLen = find(ismember(calcLengths,barsPassThresh));
% passthreshLen = barsPassThresh(calcLengths);

bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], compStr(passthreshLen), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

% compStr{passthreshLen(1)}.pos = 4169;
% super_quick_plot(1,bgAll(4),compStr(2),theoryStruct)

% allCoefs = cellfun(@(x) x.maxcoef(1),compStr);
% allLengths =300;
% 
% idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
% coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% 
% barsPassThresh = find(coefDiffs > 0);
% 
% bg = bgAll(barsPassThresh); % make gen_reference_based_assembly it work for MP
[outConsensus, coverage, pval] = gen_reference_based_assembly(bgAll(passthreshLen),compStr(passthreshLen),theoryStruct,'test11',inf);


%% Discriminative analysis (when theory bar is not known, but is known it should be e-coli or a specific genome)

goodbars = bgAll(passthreshLen);
w = 300;
[rezMax,bestBarStretch,bestLength,discSpecies] = local_alignment_assembly(goodbars, nmbp,w);