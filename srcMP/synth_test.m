% test: generate synthetic barcodes, run MP, plot

%% General parameters

SCAMP_LINE = '~/git/SCAMP/build/SCAMP';
MIN_OVERLAP_PIXELS = 300;
STRETCH_MAX = 0.05;
STRETCH_STEP = 0.01;
MIN_BARCODE_SCORE_THRESHOLD = 3;
MIN_BARGROUP_SCORE_THRESHOLD = 6;
P_COMBINE_METHOD = 'stouffer';
PRINT_TO_WINDOW = 1;

NMPX = 110;
NMPSF = 300;
NMBP = 0.22;


delete(gcp('nocreate'))
numWorkers = 24;

parpool('local',numWorkers)

%% --- Load experimental barcodes ---

TOTAL_RAND_LENGTH = 10000;
PSF_WIDTH_PIXELS = NMPSF / NMPX;
NUM_RAND_FRAGMENTS = 100;
MEAN_FRAGMENT_LENGTH = 850;
FRAGMENT_LENGTH_STD = 100;
ADDED_NOISE_MEAN = .4;
ADDED_NOISE_STD = .1;
FRAGMENT_STRETCH_STD = 0.02;
sF = [0.95:0.01:1.05];

import TempStuff.gen_random_fragments
% [barcodes, refBarcode] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD);
% if random
import Nullmodel.gen_random;
[barcodes] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,500,10)';

bitmasks = cellfun(@(x) true(size(x)), [barcodes], 'un', 0);

%% outlier code
% [~, lb, ub] = isoutlier(cell2mat(barcodes(:)'), 'mean', 'thresholdfactor', 7);
% bitmasks = cellfun(@(x) not(isoutlier(x, 'thresholdfactor', 7)) ...
%     & (x > lb & x < ub) ...
%     & ((1:length(x)) > numUncert & (1:length(x)) < length(x) - numUncert + 1), barcodes, 'un', 0);

%TODO: we should also know the correct position of these barcodes along the
%refbarcode..

barSynth = cell2struct([barcodes'; bitmasks']',{'rawBarcode','rawBitmask'},2);
lengths = cellfun(@(x) length(x),barcodes);

% barSynth = cell(1,length(barcodes));
% for i=1:length(barcodes)
%     barSynth{i}.rawBarcode = barcodes{i};
%     barSynth{i}.rawBitmask = bitmasks{i};
% end


% have to save separately..
foldSynth= 'synth';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barSynth,sF,foldSynth);


[names2, baridx2] = Core.save_long_bar_txt_skip(barSynth,1:length(barSynth),foldSynth);


%%
%%  calculates all MP and MPI
import Zeromodel.beta_ev_cdf;
NN=length(barSynth);
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
rescaleFactorPair1=cell(1,NN);
orPair1=cell(1,NN);
maxMP = zeros(1,NN);
pvalPars = zeros(2,NN);
pvals = zeros(1,NN);
tic
parfor k=1:NN
    
    % command = strcat(['~/git/SCAMP/build/SCAMP --window=' num2str(minLen) ' --input_a_file_name='...
    %    names2{k} ' --input_b_file_name=' names{k} ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
    %    ' --output_a_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
    %    ' --output_a_index_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);

    command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name='...
       names2{k} ' --input_b_file_name=' namesBar{k} ' --num_cpu_workers=28 --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);

    [status,cmdout] = system(command);

    mp1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
    mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));

    maxMP(k) = max(mp1{k});
    
    [PKS1{k},LOCS1{k},pksUnique1{k},pksUniquePos1{k},barcodePair1{k},rescaleFactorPair1{k},orPair1{k}] = ...
    Core.unique_peak_barcodes(mp1{k},mpI1{k}, 0.4,baridx2,stridx,k);

    % will also compute p-value
    a = 0.13*MIN_OVERLAP_PIXELS;
    len2 = lengths(barcodePair1{k}(1));
    len1 = lengths(k);
    n = 2*(len2-MIN_OVERLAP_PIXELS)*(len1-MIN_OVERLAP_PIXELS)/150;

    pvals(k) = 1-beta_ev_cdf(max(mp1{k}), a, 1, n, 1);
%     [p] = 

    
end
toc

mean(pvals)

%% test plot best barcode pair
[a,b] = sort(maxMP,'descend');

% import Core.plot_best_match;
ix1 = b(1);

ix = 1;
import Core.plot_best_match;
 f = plot_best_match(ix,stridx,mpI1{ix1},LOCS1{ix1},pksUniquePos1{ix1},ix1,baridx2{ix1},sF,barSynth,MIN_OVERLAP_PIXELS);


 %%
allbars = b;
[c,idPos] = sort(b);
allb = ones(1,length(allbars));
%%
bvec= [];
for ii=1:1
    if sum(allb) ==0
        break;
    end
    % now loop through first few
    btemp = allbars(find(allb));
    idx = 1;
    bar1 = btemp(idx);

    bar1
    % % now for a single barcode check the cluster barcodes
    % bar1= goodIdx(8);  % 709
    % bar1= 548;  % 709
    peaksToTry = pksUnique1{bar1};
%     allb(peaksToTry)
    % now create bargroup for bar1.
    import Core.plot_bargroup;
    [barMat{ii},bars{ii}] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barSynth,MIN_OVERLAP_PIXELS);
     %                                       ix,stridx,mpI                   ,LOCS                   ,pksUnique1                  ,pksUniquePos                   ,k,               baridx, sF, barcodeGenGood, h
    
    allb(idPos([bar1 peaksToTry])) = 0;
%     bvec{ii} = [bar1 peaksToTry];
end

%% consensus

barCons = [];
% 
% for ii=1
%     
% end

stMat = min(cellfun(@(x) min(x{1}),barMat{ii}));
stopM = max(cellfun(@(x) max(x{1}),barMat{ii}));

consensusMat = nan(length(barMat{ii}),stopM-stMat+1);
% 
for i=1:length(barMat{ii})
    consensusMat(i,barMat{ii}{i}{1}-stMat+1:barMat{ii}{i}{1}-stMat+length(barMat{ii}{i}{2})) = barMat{ii}{i}{2};
    
end

figure,plot(consensusMat')
% 
% 
% barMat{1}{1}{1}
% 
