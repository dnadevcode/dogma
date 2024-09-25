% Script to check if we match the experimental barcode to correct sequence,
% based on which sample does the experimental barcode come from (and the
% sample is a priori sequenced, so sequence is known)
ffile = {'/export/scratch/albertas/data_temp/bargrouping/ecoli/ecoli_2_old_1/',...
    '/export/scratch/albertas/data_temp/bargrouping/ecoli/ecoli_2_old_2/',...
    '/export/scratch/albertas/data_temp/bargrouping_selected/ecoli_1/sample_2/'};

ix=2
import Core.load_chrom_data;
[bgAll, bG, kymoStructs] = load_chrom_data(ffile{ix});


numWorkers = 30; % num workers, in one node there are 30
% minLen = [150:50:3000]; % min to max % could take more points for more accurate..
sets.comparisonMethod = 'mass_pcc';
sF = 0.9:0.01:1.1;
nmPerPx = 110;

%todo: rename:  018_final_polish - Sample2_EF365_DA70365
%               contigs_EF544 - Sample5_EF544_DA?
%               
fastas = {'018_final_polish.fasta','DA32087.fasta'};%'contigs_EF544.fasta'};

import Thry.gen_theoretical;
fastaFiles = 1:length(fastas);

% if for all nm/bp values
nmbpvals = 0.27; % should use mpMax to get this "fully" correct
[theoryStructNew,~,barcodeGenT] = gen_theoretical(fastas(fastaFiles),nmbpvals,0,nmPerPx); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly

% load('/export/scratch/albertas/data_temp/bargrouping/ecoli/ecoli_2_new/session_data_3.mat');

barGenRun = bgAll;
% barGenRun = barcodeGen;

w = [];
[rezMax,bestBarStretch,bestLength,rezOut] = local_alignment_assembly(theoryStructNew, barGenRun,w);


super_quick_plot(1,barGenRun,rezOut{1},theoryStructNew)
% rezOut{1}.rezMax{2}{3}


import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name(theoryStructNew(1:2));
% 
idx = 1;
import Core.discrim_true_positives;
[truePositives,discSpecies,discAll,allSpecies,refNums,signMatch] =...
    discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);



% figure,histogram(cell2mat(refNums))
% refs = cell2mat(refNums)
refBest = cellfun(@(x) x(1),refNums);
[sum(refBest==1)/length(refBest) sum(refBest==2)/length(refBest)]


%%
%% compare and get discriminative stuff
nmvals = 0.22:0.005:0.3;
vals = [];
for nmbpvals = nmvals
    
    % if for all nm/bp values
    % nmbpvals = 0.225; % should use mpMax to get this "fully" correct
    [theoryStructNew,~,barcodeGenT] = gen_theoretical(fastas(fastaFiles),nmbpvals,0,nmPerPx); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly
    
    
    barGenRun = bgAll;
%     w = 300;
    [rezMax,bestBarStretch,bestLength,rezOut] = local_alignment_assembly(theoryStructNew, barGenRun,[]);

    [truePositives,discSpecies,discAll,allSpecies,refNums,signMatch] =...
        discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);

    
    refBest = cellfun(@(x) x(1),refNums);
    vals(end+1,:) = [sum(refBest==1)/length(refBest) sum(refBest==2)/length(refBest)];
end

figure,plot(nmvals,vals)
legend({'Ecoli1','Ecoli2'})