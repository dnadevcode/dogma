timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% input parameters
pthresh = 0.001; % small threshold for creating bargroup
minAv = 0; % minimum number of barcodes in consensus


NUM_RAND_FRAGMENTS = 100;% number of random fragments
PSF_WIDTH_PIXELS = 4.7; % psf
MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
FRAGMENT_LENGTH_STD = 100; %std of length of fragment
ADDED_NOISE_MEAN = 0.0; % additive noise mean
ADDED_NOISE_STD = 0.05; % additive noise std
FRAGMENT_STRETCH_STD = 0;%0.02; % fragment length-rescale factor std 0.02
IS_CIRCULAR = 1; % whether circular barcode
minL = 150; % minimum length

%% generate synthetic data
import Nullmodel.gen_rand;
[barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr] = gen_rand(...
    NUM_RAND_FRAGMENTS, PSF_WIDTH_PIXELS,MEAN_FRAGMENT_LENGTH,FRAGMENT_LENGTH_STD,...
    ADDED_NOISE_MEAN,ADDED_NOISE_STD,FRAGMENT_STRETCH_STD,IS_CIRCULAR,minL);

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);

%% theory self-similarity

%% visualize true output
import Plot.plot_best_pos;
% fig1 = plot_best_pos([], synthStr, [], [], [],theoryStruct{1}.length,1);%1/bpPx*10^6

%% visualize best aligment
idxPair = 1;
import Core.plot_match_pcc;
[sortedVals, sortedIds,pscores] = sorted_scores(synthStr2);
[curSink,curSource] = ind2sub(size(synthStr2),sortedIds(idxPair));
[f] = plot_match_pcc(barStruct, synthStr2,curSink,curSource,barStruct);

%% create barcoding islands
import Core.create_barset;
[barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = create_barset(pscores,sum(~isnan(pscores(:))), barcodeGen, synthStr2,synthStr);

%% plot barcode islands
nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
import Plot.islandsPosStruct;
posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
import Plot.plot_best_pos;
% fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], synthStr(cell2mat(barcodeIslands)), [], [], [],10000,1);%1/bpPx*10^6


%% reference based assembly: TODO: possible improvements?
[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,synthStr,theoryStruct,timestamp);
% [outConsensus, coverage] = gen_assembly(barcodeGen,synthStr,cumsum(repmat(10000,1,1)),timestamp);

% assembly w/o reference / bars are shuffled..
% [outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(barcodeGen(cellfun(@(x) x.barid,posStruct)),posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp);

%% plot bargrouping generated island:
% pthresh = 0.00001;
import Plot.plot_island;
% import Core.simple_bargroup;
% plot_island(posStruct,islandElts, barcodeGen(cellfun(@(x) x.barid,posStruct)),1);


%% Up until now was ground truth, i.e. the data that we have. 
% Now do the calculations. 1) calculate scores:
sF = 1;%0.95:0.01:1.05;
minOverlap = 150; % test for best min-overlap
scorethresh = 0.5; % instead use Stouffer like for HC
tic
import Core.calc_overlap_pcc_sort_m;
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap,1);
toc

% idxPair = 7;
% import Core.plot_match_pcc;
% [sortedValsC, sortedIdsC,~] = sorted_scores(overlapStruct);
% [curSink,curSource] = ind2sub(size(overlapStruct),sortedIdsC(idxPair));
% [f] = plot_match_pcc(barStruct, overlapStruct,curSink,curSource,barStruct);

    
res.overlapStruct = overlapStruct;

% timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
mkdir(strcat('output',timestamp));

% to check local similarity/remove false positives. Skip for now for speed
% [overlapStructMP] = calc_overlap_mp(barcodeGen',sF, minOverlap-50,timestamp)

%% now barset      
import Zeromodel.pval_for_struct;
pscores = pval_for_struct(overlapStruct,scorethresh,0.05);

pthresh = 0.001;
pscores(pscores>pthresh) = nan;

psc = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1
filtM = triu(psc);
filtM(filtM==0) = nan;

    
import Core.create_barset;
[barcodeIslands, barcodeIslandsData, badData, badDataInfo, barIslands] = create_barset(-pscores,sum(~isnan(filtM(:))),barcodeGen,overlapStruct);
nonEmpty = find(cellfun(@(x) ~isempty(x),barcodeIslands));
import Plot.islandsPosStruct;
posStructE = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));

% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], posStructE, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6


    % assembly w/o reference / bars are shuffled..
[outConsensus, coverage, consensus,islandEltsE, islandPx] = gen_assembly(barcodeGen(cellfun(@(x) x.barid,posStructE)),posStructE,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp,minAv);
import Plot.plot_island;

[outConsensus] = plot_island(posStructE,islandEltsE, barcodeGen(cellfun(@(x) x.barid,posStructE)),1);
% plot_island(posStruct,islandElts, barcodeGen,1);

% TP / FN rate
% calc_discovery_rates(outConsensus,theoryStruct,synthStr)
[posDif] = calc_discovery_rates(outConsensus,theoryStruct,synthStr,posStructE,islandEltsE);

% 
% function [sortedVals, sortedIds,pscores] = sorted_scores(synthStr2)
% 
%     psc = reshape([synthStr2(:).score], size(synthStr2,1),size(synthStr2,2)); % from bg_test_1
%     filtM = triu(psc);
%     filtM(filtM==0) = nan;
%     pscores = filtM(:);
%     [sortedVals, sortedIds] = sort(pscores,'desc','MissingPlacement','last');
% 
% end
