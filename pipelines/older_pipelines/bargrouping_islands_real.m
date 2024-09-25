timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%%
%%
matFileLoc = 'C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_09_27_22\barcodeGen993.mat';
fastaFileLoc = 'C:\Users\Lenovo\postdoc\DATA\Chromosome\YEAST\*.fasta';

%%
minLen = 200; % minimum length
numF = 10; % number of frames 
sF = 1;%0.95:0.01:1.05;
minOverlap = 250;
scorethresh = 0.7; % instead use Stouffer like for HC

    
matFile = 1;
if matFile
    load(matFileLoc);
else
    testSet = {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/'};
    % testSet ={ '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
    %     '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
    synth = 0; % if synthetic

    % minOverlap = 150;
    import Zeromodel.prep_data;
    [barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

    barcodeGen = barcodeGen1;
end

% NUM_RAND_FRAGMENTS = 100;% number of random fragments
% PSF_WIDTH_PIXELS = 4.7; % psf
% MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
% FRAGMENT_LENGTH_STD = 100; %std of length of fragment
% ADDED_NOISE_MEAN = 0.4; % additive noise mean
% ADDED_NOISE_STD = 0.1; % additive noise std
% FRAGMENT_STRETCH_STD = 0.02;%0.02; % fragment length-rescale factor std 0.02
% IS_CIRCULAR = 1; % whether circular barcode
% minL = 150;
% 
% %% generate synthetic data
% import Nullmodel.gen_rand;
% [barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr] = gen_rand(...
%     NUM_RAND_FRAGMENTS, PSF_WIDTH_PIXELS,MEAN_FRAGMENT_LENGTH,FRAGMENT_LENGTH_STD,...
%     ADDED_NOISE_MEAN,ADDED_NOISE_STD,FRAGMENT_STRETCH_STD,IS_CIRCULAR,minL);

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

%% Check for barcodes which are periodic: run self similarity
import Core.find_periodic_barcodes;
[valsPassingThresh,pccValsMax,mp1,mpI1] = find_periodic_barcodes(barStruct,timestamp)

barStruct = barStruct(find(valsPassingThresh));

%% check and merge similar neighbours (in case this was not identified via DBM)

%% Mask 'bad' regions

barcodes = arrayfun(@(x) barStruct(x).rawBarcode,1:length(barStruct),'un',false);
bitmasks = cellfun(@(x, y) y & abs(x - nanmean(x(y))) <= 2.5 * nanstd(x(y)),...
    arrayfun(@(x) barStruct(x).rawBarcode,1:length(barStruct),'un',false),...
    arrayfun(@(x) barStruct(x).rawBitmask,1:length(barStruct),'un',false), 'un', 0);

%% visualize true output
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], synthStr, [], [], [],theoryStruct{1}.length,1);%1/bpPx*10^6

%% visualize best aligment
% idxPair= 1;
% import Core.plot_match_pcc;
% [sortedVals, sortedIds,pscores] = sorted_scores(synthStr2);
% [curSink,curSource] = ind2sub(size(synthStr2),sortedIds(idxPair));
% [f] = plot_match_pcc(barStruct, synthStr2,curSink,curSource,barStruct);

%% create barcoding islands
% import Core.create_barset;
% [barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = create_barset(pscores,sum(~isnan(pscores(:))), barcodeGen, synthStr2,synthStr);

%% plot barcode islands
% nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
% import Plot.islandsPosStruct;
% posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], synthStr(cell2mat(barcodeIslands)), [], [], [],10000,1);%1/bpPx*10^6


%% comp to thry 
fastaFile =[]; % todo:test

fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};


nmbp = 0.216;
nmpx = 110;
psffac = 1;
nmbp = nmbp*psffac;
nmpx = nmpx*psffac;

import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmpx);

sets.comparisonMethod = 'mass_pcc';

sF = 0.9:0.01:1.1;
[comparisonStruct] = compare_to_t(barcodeGen,theoryStruct,sF,sets)


bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

% MP:
[compStr] = compare_to_t_mp(barcodeGen,theoryStruct,1,150)

%% reference based assembly: TODO: possible improvements?: for this: compare against thr
[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,comparisonStruct,theoryStruct,timestamp);
[outConsensus, coverage] = gen_assembly(barcodeGen,comparisonStruct,cumsum(repmat(10000,1,1)),timestamp);

% assembly w/o reference / bars are shuffled..
[outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(barcodeGen(cellfun(@(x) x.barid,posStruct)),posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp);

%% plot bargrouping generated island:
% pthresh = 0.00001;
% import Plot.plot_island;
% import Core.simple_bargroup;
% plot_island(posStruct,islandElts, barcodeGen(cellfun(@(x) x.barid,posStruct)),1);

%% calculate scores:
sF = 1;
tic
import Core.calc_overlap_pcc_sort_m;
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap,1);
toc
    
res.overlapStruct = overlapStruct;

% timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
mkdir(strcat('output',timestamp));

% to check local similarity/remove false positives. Skip for now for speed
% [overlapStructMP] = calc_overlap_mp(barcodeGen',sF, minOverlap-50,timestamp)

%% now barset      
import Zeromodel.pval_for_struct;
pscores = pval_for_struct(overlapStruct,scorethresh,0.05);


psc = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1
filtM = triu(psc);
filtM(filtM==0) = nan;

    
import Core.create_barset;
[barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = create_barset(pscores,sum(~isnan(filtM(:))),barcodeGen',overlapStruct);

    % assembly w/o reference / bars are shuffled..
nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
import Plot.islandsPosStruct;
posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
import Plot.plot_best_pos;
fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6


    [outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(barcodeGen(cellfun(@(x) x.barid,posStruct)),posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp);
import Plot.plot_island;

plot_island(posStruct,islandElts, barcodeGen(cellfun(@(x) x.barid,posStruct)),1);
