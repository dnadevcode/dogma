addpath(genpath('/home/avesta/albertas/reps/lldev'))
addpath(genpath('/home/avesta/albertas/reps/hca'))
addpath(genpath('/home/avesta/albertas/reps/discriminative_sequences'))

addpath('/export/scratch/albertas/data_temp/bargrouping/ecoli/FASTAS/')
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/'))
addpath(genpath('/home/avesta/albertas/reps/bargroupingprototype'))

% GUI can be loaded via bargi_gui


% EF544 is DA68355. EF365 is DA68332.

%% All figs:
%
% Main text:
% Fig1 - overlap + evaluationcco
% Fig2 - contig
% Fig3/Table 1 - coverage from each day experiment
% Fig4 - overlap against theory

% Supplementary figs    
% FigS1 - local similarities threshold
% FigS2 - success rate wrt nm/bp ratio for individual day experiment
% FigS10 - local similarities threshold validation on synthetic data -
%so validation for S2
% FigS3 - distance between two overlapping barcodes vs. corresponding
% distance of reference mapped barcodes. 
%

%% data loading:

% synthetic data (mostly for evaluation)
import Core.load_synth_data;
[bgAll, bG, synthStr, synthStr2, theoryStruct,setsGen] = load_synth_data();


dataFold = {'/export/scratch/albertas/data_temp/bargrouping_selected/'};
dataSets = {'ecoli_S2/','ecoli_S5/','ecoli_dEC','synth','ecoli_test/','ecoli_test2/','ecoli_test3/','ecoli_EF365_large/'}; % different available datasets
fastas = {'018_final_polish.fasta','DA32087.fasta','DA68335.fasta'};

% 2-nd probably not 1
theoryIdxs = {1, 3, 2, nan, 2 , 2, 3,1}; % known theory indexes

dataSetIdx = 3; % 5,... test datasets (part of the data)
import Core.load_chrom_data;
[bgAll, bG, kymoStructs] = load_chrom_data(fullfile(dataFold{1},dataSets{dataSetIdx}),'kymo');

 
kymoParams.nmbp = cellfun(@(x) x.nmbp,bgAll);
kymoParams.nmpx = cellfun(@(x) x.nmpx,bgAll);
import TempStuff.rescale_barcodes;
barcodeGen2 = rescale_barcodes(bgAll, kymoParams, 110, 0.2);


% If needs to be generated:
% userDir = '/export/scratch/albertas/data_temp/bargrouping/ecoli/ecoli_2_new/1219/';
% [bgAll,bG,kymoStructs] = load_chrom_data(userDir, 'tiff', 0, 20);


% Evaluation for each data-set separately
fig_evaluation_individual_days

% evaluation for all together
fig_evaluation

%
fig_evaluation_simulated


% contigs
fig_contigs
fig_contigs_combined
%
fig_contigs_simulated

%% Validation
figs_validation_on_synthetic_data % figure S10.eps
figs_validation_on_real_data % figure S10.eps

fig_window_width


% 
% figs_test_on_real_data % here we use data from experiment we don't have corresponding barcode from
  

%%
% Figure 1
fig_example_1

% For Fig3/Fig4
example_plot

%% All paper figures (manuscript version 6)

% Figure 1
fig1_ex

% Figure 2
% assembly_schematics.drawio

final_figure_3


%% Supl Fig S1 (FigS1.eps)
import Nullmodel.PvalScripts.pval_test_loop_pcc;
pval_test_loop_pcc;
%% Supl Fig S2 (FigS2.eps)

import Nullmodel.PvalScripts.pval_test_synth_local;
pval_test_synth_local

%% Supl Fig S3
import Validation.length_Validation;
length_Validation