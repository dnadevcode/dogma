% this script is chromosome assembly script vs reference, where both input
% and reference are synthetic barcodes.

% add hca to path for theory calculation and pcc comparison
% addpath(genpath('C:\Users\Lenovo\git\hca'));
addpath(genpath('/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/hca'))
addpath(genpath(pwd));

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

% rootdir = 'C:\Users\Lenovo\postdoc\DATA\Frag\all2';

%% gen synthetic
import Nullmodel.gen_random_fragments;

% parameters
TOTAL_RAND_LENGTH = 10000;
PSF_WIDTH_PIXELS = 2.7;
NUM_RAND_FRAGMENTS = 100;
MEAN_FRAGMENT_LENGTH = 850;
FRAGMENT_LENGTH_STD = 100;
ADDED_NOISE_MEAN = 0.4;
ADDED_NOISE_STD = .1;
FRAGMENT_STRETCH_STD = 0.02;
sF = [0.9:0.01:1.1];
IS_CIRCULAR=1;
    
% generate barcodes
import Nullmodel.gen_random_fragments;
[barcodeGen,refBarcode,origPos, origFlip,origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD,IS_CIRCULAR);

% barSynth = cell2struct([barcodes'; bitmasks']',{'rawBarcode','rawBitmask'},2);
barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);

% lengths of barcode masks
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);

% plot without calculating scores. or can be used to pass after calculating
% using only good
import Nullmodel.plot_random_fragments
[outConsensus] = plot_random_fragments(barcodeGen,1:length(barcodeGen),refBarcode,origPos,1./origStr,origFlip'+1);


[a,msg] = mkdir(strcat('output',timestamp));
i=1;theoryStruct=[];
theoryStruct{i}.filename = fullfile(strcat('output',timestamp),'theory_barcode.txt');
fileID = fopen(theoryStruct{i}.filename,'w');
fprintf(fileID,strcat([' %2.' num2str(14) 'f ']), refBarcode);
fclose(fileID);
theoryStruct{i}.meanBpExt_nm = 0.3;
theoryStruct{i}.psfSigmaWidth_nm = 300;
theoryStruct{i}.length = length(refBarcode);
theoryStruct{i}.isLinearTF = 0;
theoryStruct{i}.name = 'Synthetic theory';


% todo: compare in the same way as MP. 
% can create one long file from theoryStruct.
rezMax=[];bestBarStretch=[];bestLength=[];
for i=1:length(theoryStruct)
 tic
 import CBT.Hca.Core.Comparison.on_compare;
[rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGen,theoryStruct{i},'mass_pcc',sF,[],50);
toc
end

comparisonStructAll = rezMax;
for i=1:length(comparisonStructAll)
    for j=1:length(bestBarStretch{i})
        comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
        comparisonStructAll{i}{j}.length = bestLength{i}(j);
    end
end
import CBT.Hca.Core.Comparison.combine_theory_results;
[comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);

%% create consensus
[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,comparisonStruct,theoryStruct);

%% consensus specific maps
%
barIdx = [4,8];
[outConsensus] = gen_reference_based_assembly(barcodeGen(barIdx),comparisonStruct(barIdx),theoryStruct);

% outConsensus2=circshift(outConsensus,[0,10]);
% f=figure,plot(find(~isnan(nanmean(outConsensus2(2:end,:)))),outConsensus2(2:end,~isnan(nanmean(outConsensus2(2:end,:))))'); 

%% now de-novo assembly (in the same script)

delete(gcp('nocreate'))
% numWorkers = 28;
% 
% parpool('local',numWorkers)
%%
c = parcluster('kebnekaise');
c.AdditionalProperties.AccountName = 'snic2021-5-132';
c.AdditionalProperties.WallTime = '72:00:00';

numWorkers = 28;
% addpath(genpath('C:\Users\Lenovo\git\bargroupingprototype'))
% cd C:\Users\Lenovo\git\bargroupingprototype


% SCAMP_LINE = '.\..\SCAMP\build\Release\SCAMP.exe'; % windows
SCAMP_LINE = '~/git/SCAMP/build/SCAMP'; % how to call scamp from bash

MIN_OVERLAP_PIXELS = 300;
STRETCH_MAX = 0.0;
STRETCH_STEP = 0.01;
MIN_BARCODE_SCORE_THRESHOLD = 3;
MIN_BARGROUP_SCORE_THRESHOLD = 6;
P_COMBINE_METHOD = 'stouffer';
PRINT_TO_WINDOW = 1;

NMPX = 110;
NMPSF = 200;
NMBP = 0.22;
sF = 0.9:0.01:1.1;

%%

% barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
%     cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);
% have to save separately..
foldSynth= 'barcoli';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);

% parpool(c,10)

 parpool('local',numWorkers)
[names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);
% delete(gcp('nocreate'))




% parpool(c,numWorkers)

% batchPoolSize = 12;

% c2 = parcluster('local');
% delete(gcp('nocreate'))
% parpool(c2)
%%  calculates all MP and MPI
tic
import Core.calc_overlap;
[mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,names2,namesBar);
toc
delete(gcp('nocreate'))

% tic
% j = c2.batch(@Core.calc_overlap, 3, {barStruct,timestamp,SCAMP_LINE,MIN_OVERLAP_PIXELS,numWorkers,names2,namesBar,c2}, 'Pool', batchPoolSize);
% j.wait;
% scoreRegistry = j.fetchOutputs{:};
% toc
% 
% j = c.batch(@Core.calc_overlap, 3, {barStruct,timestamp,SCAMP_LINE,numWorkers,names2,namesBar}, 'Pool', batchPoolSize);
% 
% %%
% 
% j.wait;
% scoreRegistry = j.fetchOutputs{:};
% toc
% parpool('local',28)
NN=length(mpI1);
PKS1=cell(1,NN);
LOCS1=cell(1,NN);
pksUnique1 = cell(1,NN);
pksUniquePos1=cell(1,NN);
barcodePair1=cell(1,NN);
rescaleFactorPair1=cell(1,NN);
orPair1=cell(1,NN);
pvals =zeros(1,NN);
for k=1:NN
    k
    
    [PKS1{k},LOCS1{k},pksUnique1{k},pksUniquePos1{k},barcodePair1{k},rescaleFactorPair1{k},orPair1{k}] = ...
    Core.unique_peak_barcodes(mp1{k},mpI1{k}, 0.3,baridx2,stridx,k);
% 
% %     % will also compute p-value
%     a = 0.122*MIN_OVERLAP_PIXELS;
%     len2 = lengths(barcodePair1{k}(1));
%     len1 = lengths(k);
%     n = 2*(len2-MIN_OVERLAP_PIXELS)*(len1-MIN_OVERLAP_PIXELS)/100;
% % 
%     pvals(k) = 1-Zeromodel.beta_ev_cdf(maxMP(k), a, 1, n, 1);
% % %

    
end

% calculate full overlap PCC scores
import Core.mp_to_full_overlap;
[PCC_OVERLAP,len1,len2] = mp_to_full_overlap(stridx,pksUnique1,pksUniquePos1,LOCS1,mpI1,baridx2, MIN_OVERLAP_PIXELS, barStruct,sF);

%%
%%
[a,b] = sort(maxMP,'descend');

% import Core.plot_best_match;
ix1 = b(1)
% % ix1= 3
% ix1=195;
ix = 1;
pos = LOCS1{ix1}(pksUniquePos1{ix1}(ix)); % position on MP.

import Core.plot_best_match;
[f,~,~,~,~,~,bar1Name] = plot_best_match(ix,stridx,mpI1{ix1},LOCS1{ix1},pksUniquePos1{ix1},ix1,baridx2{ix1},sF,barStruct,MIN_OVERLAP_PIXELS);



barIdx = [ix1  bar1Name];

% if theory generated..
[outConsensus] = gen_reference_based_assembly(barcodeGen(barIdx),comparisonStruct(barIdx),theoryStruct);

%%
import Core.overlap_graph
G = overlap_graph

import Core.layout_consensus
layout_consensus
