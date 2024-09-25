% validation of individual overlaps for barcodes. Depends if they match
% correctly on the theory..


%% Generate theory
% First generate theory:
sF = 0.7:0.025:1.3;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

fastas = {'018_final_polish.fasta','DA32087.fasta','DA68335.fasta'};
theoryIdxs = {1, 3, 2, nan, 2 , 2, 3, 1}; % known theory indexes
dataSetIdx = 3; % 5,... test datasets (part of the data)

nmPerPx = 110;
nmbp = 0.2; % ? 

% quickly calculate theory for given single fasta file
[theoryStructRev,theoryStruct,bT] = prep_thry_for_local_comp(fastas(theoryIdxs{dataSetIdx}), nmbp, nmPerPx, 1);


%% Load bars
load('2all_2023-11-28_09_43_54.mat')
% load('3all_2023-11-24_12_38_55.mat');
pval =  0.01;
alphaNu = 0.08;
sets.minOverlap = 300;
% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pval); %localCCStruct

%%1) Run comparison (all barcodes). First convert to same nm/bp ratio
%% (assuming nm/bp from lambda's is correct)
% bars = barcodeGen(barIslands{idx}); 
% barx = barIslands{idx};
sets.comparisonMethod = 'mass_pcc';
sets.nuF = 0.08;
[compI,rezI,~] = compare_to_t(barcodeGen,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?

[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,compI,{theoryStruct},'test11',inf);



N = 1000;
import Validation.calc_distance_scaled; % need to properly scale things when calculating distances for validation
[posShift,sinkSource] = calc_distance_scaled(oS,sortedIds(1:N), compI,theoryStruct.length, N);
figure,histogram(posShift<50)
%  
figure,plot(posShift<50)

idx = find(posShift<50);

allBadMatch = sinkSource(idx,:);


%% plot individual overlap vs theory

idxPair = 1;
 

[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
% 
% lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
lenThry = theoryStruct.length;
import Plot.pair_evaluation_plot;
[f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,compI([curSink curSource]),lenThry);


import Core.comparisonsres_to_struct; 
[rezStruct2] = comparisonsres_to_struct(barcodeGen, compI,lenThry);


import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
[f] = pair_evaluation_with_ground_truth_plot(barcodeGen, oS,curSink,curSource,rezStruct2,lenThry);


sets.comparisonMethod = 'mpnan';
sets.w = 300;
sets.nuF = 0.08;
[compI2,rezI2,~] = compare_to_t(barcodeGen([curSink,curSource]),theoryStruct,sF,sets); % maybe calculate Stouffer score here also?
[compI2] = Core.mp_to_pcc_pos(compI2); % convert to pcc plotable

compI2{1}
compI2{2}


% theoryStruct.filename = '1';
import CBT.Hca.UI.Helper.plot_any_bar;
plot_any_bar(1,barcodeGen(curSink),{{compI{curSink}}},theoryStruct,1);


import CBT.Hca.UI.Helper.plot_any_bar;
plot_any_bar(1,barcodeGen(curSource),{{compI{curSource}}},theoryStruct,1);




ix = 1;
import CBT.Hca.UI.Helper.plot_any_bar;
plot_any_bar(1,barcodeGen(curSink),{{compI2{1}}},theoryStruct,1);


% ix = 2;
import CBT.Hca.UI.Helper.plot_any_bar;
plot_any_bar(1,barcodeGen(curSource),{{compI2{2}}},theoryStruct,1);


% mp to regular?

