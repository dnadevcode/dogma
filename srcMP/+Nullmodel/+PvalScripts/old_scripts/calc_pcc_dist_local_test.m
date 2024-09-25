% calc_pcc_dist_local_test (previously pval_test)

%{
 This just calculates PCCs for given lengths, so not the full pvalue test

 P-value depends on these parameters: 
- overlap length
- number of attempts, re-scaling factors
- psf.
if re-scaling factors and psf have minor relevance , then it would
 be mainly the other two parameters that we would need to tune.


% for the analysis, the total number of fitting attempts will be important,
% not the lengths of individual barcodes. For the analysis we take common
% length, and short barcode
%}



import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 1000;
minOverlapPx = 300;
PSF_WIDTH_PIXELS = 300/110;
randLenPx = 500;
randLengthGapPx = 500;

sF = 0.95:0.01:1.05;

[~ , barSynth] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,randLenPx,0);
%
[~, barSynth2] = gen_random(1,PSF_WIDTH_PIXELS,randLenPx,randLengthGapPx);



% have to save separately..
foldSynth= 'synthtest';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barSynth2,sF,foldSynth);

% now save the same barcodes as long vectors since we will compare against
% these.\todo: possibiility of using one long barcode?
[names2, baridx2] = Core.save_long_bar_txt_skip(barSynth,inf,foldSynth);



%%  calculates all MP and MPI
NN=length(names2);
out='output';
mkdir(out);
mp1 = cell(1,NN);
mpI1 = cell(1,NN);

maxPCC=cell(1,NN);
tic
for k=1:NN
    maxPCC{k} = nan(1,NN);
    command = strcat(['SCAMP' ' --window=' num2str(minOverlapPx) ' --input_a_file_name='...
       names2{k} ' --input_b_file_name=' namesBar{k} ' --num_cpu_workers=28 --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);

    [status,cmdout] = system(command);

    mp1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));

    [PKS,LOCS] = sort(mp1{k},'descend');
    LOCS = LOCS(~isnan(PKS));
    PKS = PKS(~isnan(PKS));
    barNum = baridx2{k}(LOCS);
    [LOCSUIQUE,LOCSUNIQUEPos] = unique(barNum,'stable');
    maxPCC{k}(LOCSUIQUE) = PKS(LOCSUNIQUEPos);    
end
toc


cellfun(@(x) nanmean(x),maxPCC)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Old code

%{

% test: generate p-values for bargroups
SCAMP_LINE = '~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
delete(gcp('nocreate'))
numWorkers = 24;

parpool('local',numWorkers)


% overlap length in general will be fixed for the analysis.
minOverlapPx = 300;


% for the analysis, the total number of fitting attempts will be important,
% not the lengths of individual barcodes. For the analysis we take common
% length, and short barcode
import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 1000;
PSF_WIDTH_PIXELS = 300/110;
randLenPx = 500;
randLengthGapPx = 500;

sF = 0.95:0.01:1.05;

[barcodes] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,randLenPx,0);
%
bitmasks = cellfun(@(x) true(size(x)), [barcodes], 'un', 0);


% create bar structure 
barSynth = cell2struct([barcodes; bitmasks]',{'rawBarcode','rawBitmask'},2);


[barcodes2] = gen_random(1,PSF_WIDTH_PIXELS,randLenPx,randLengthGapPx);
bitmasks2 = cellfun(@(x) true(size(x)), [barcodes2], 'un', 0);

barSynth2 = cell2struct([barcodes2; bitmasks2]',{'rawBarcode','rawBitmask'},2);


% barSynth = [cell2struct(barcodes,{'rawBarcode'},1);cell2struct(bitmasks,{'bitmasks'},1)];


% have to save separately..
foldSynth= 'synthtest';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barSynth2,sF,foldSynth);

% now save the same barcodes as long vectors since we will compare against
% these.\todo: possibiility of using one long barcode?


[names2, baridx2] = Core.save_long_bar_txt_skip(barSynth,inf,foldSynth);


% barSynth = cellfun(
% barSynth = cell(1,length(barcodes));
% for i=1:length(barcodes)
%     barSynth{i}.rawBarcode = barcodes{i};
% end

%%
%%  calculates all MP and MPI
NN=length(names2);
out='output';
mkdir(out);
% parpool('local',28)
mp1 = cell(1,NN);
mp2 = cell(1,NN);
PKS1=cell(1,NN);
LOCS1=cell(1,NN);
pksUnique1 = cell(1,NN);
pksUniquePos1=cell(1,NN);
barcodePair1=cell(1,NN);
maxPCC=cell(1,NN);
orPair1=cell(1,NN);
maxMP = zeros(1,NN);
tic
parfor k=1:NN
    maxPCC{k} = nan(1,NN);
    command = strcat([SCAMP_LINE ' --window=' num2str(minOverlapPx) ' --input_a_file_name='...
       names2{k} ' --input_b_file_name=' namesBar{k} ' --num_cpu_workers=28 --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);

    [status,cmdout] = system(command);

    mp1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));

    
%     baridx2{k}
    [PKS,LOCS] = sort(mp1{k},'descend');
    LOCS = LOCS(~isnan(PKS));
    PKS = PKS(~isnan(PKS));
    barNum = baridx2{k}(LOCS);
    [LOCSUIQUE,LOCSUNIQUEPos] = unique(barNum,'stable');
    maxPCC{k}(LOCSUIQUE) = PKS(LOCSUNIQUEPos);

% 
%     
%     maxMP(k) = max(mp1{k});
%     
%     % now we find max against each barcode
%     
%     [PKS1{k},LOCS1{k},pksUnique1{k},pksUniquePos1{k},barcodePair1{k},rescaleFactorPair1{k},orPair1{k}] = ...
%     Core.unique_peak_barcodes(mp1{k},mpI1{k}, 0.7,baridx2,stridx,k);

    
end
toc

cellfun(@(x) nanmean(x),maxPCC)

cell2mat(maxPCC)


pccMat = cell2mat(maxPCC');


% import TempStuff.gen_random_fragments
% [barcodes, refBarcode] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD);

%}
