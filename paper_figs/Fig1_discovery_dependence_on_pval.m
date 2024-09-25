timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%%
accuracyPx = 5;
pthresh = 0.0000001; % for creating bargroup

minAv = 0; % minimum number of barcodes in consensus
sF = 1;%0.95:0.01:1.05;
minOverlap = 150; % test for best min-overlap
scorethresh = 0.5; % thresh for pcc, irrelevant if using small pthresh

pvalparNu = 0.05;   % pval par. depends on PSF (calculate relation separately)


NUM_RAND_FRAGMENTS = 50;% number of random fragments
PSF_WIDTH_PIXELS = 4.7; % psf
MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
FRAGMENT_LENGTH_STD = 100; %std of length of fragment
ADDED_NOISE_MEAN = 0.4; % additive noise mean
ADDED_NOISE_STD = 0.1; % additive noise std
FRAGMENT_STRETCH_STD = 0;%0.02; % fragment length-rescale factor std 0.02
IS_CIRCULAR = 1; % whether circular barcode
minL = 150; % minimum length / same as minoverlap?


% function [outputArg1,outputArg2] = Fig1_discovery_dependence_on_pval(inputArg1,inputArg2)




%% generate synthetic data
import Nullmodel.gen_rand;
[barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr] = gen_rand(...
    NUM_RAND_FRAGMENTS, PSF_WIDTH_PIXELS,MEAN_FRAGMENT_LENGTH,FRAGMENT_LENGTH_STD,...
    ADDED_NOISE_MEAN,ADDED_NOISE_STD,FRAGMENT_STRETCH_STD,IS_CIRCULAR,minL);

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);


%% calc

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
pscores = pval_for_struct(overlapStruct,0,pvalparNu);

plist = 0.0000001:0.01:0.2;
tpr = zeros(1,length(plist));
tpr2 = zeros(1,length(plist));

pre = zeros(1,length(plist));

accuracyPx = 1;

% totNpositives

for i = 1:length(plist); % for creating bargroup
    pthresh = plist(i);
    pscr = pscores;
    pscr(pscr>pthresh) = nan;

    psc = reshape([pscr], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1
    filtM = triu(psc);
    filtM(filtM==0) = nan;


    import Core.create_barset;
    [barcodeIslands, barcodeIslandsData, badData, badDataInfo, barIslands] = create_barset(-filtM,sum(~isnan(filtM(:))),barcodeGen,overlapStruct);
    nonEmpty = find(cellfun(@(x) ~isempty(x),barcodeIslands));
    import Plot.islandsPosStruct;
    posStructE = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));

    % import Plot.plot_best_pos;
%     fig1 = plot_best_pos([], posStructE, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6


        % assembly w/o reference / bars are shuffled..
    [outConsensus, coverage, consensus,islandEltsE, islandPx] = gen_assembly(barcodeGen(cellfun(@(x) x.barid,posStructE)),posStructE,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp,minAv);
    import Plot.plot_island;

%     [outConsensus] = plot_island(posStructE,islandEltsE, barcodeGen(cellfun(@(x) x.barid,posStructE)),1);
%     plot_island(posStruct,islandElts, barcodeGen,1);

    % TP / FN rate
    % calc_discovery_rates(outConsensus,theoryStruct,synthStr)
    [posDif] = calc_discovery_rates(outConsensus,theoryStruct,synthStr,posStructE,islandEltsE);

    allvals = cell2mat(posDif');
    circCorrected  = min(abs(allvals),abs(abs(allvals)-length(refBarcode)));
    TP = sum(circCorrected < accuracyPx);
    FN = length(allvals) - TP;
    FP =  sum(circCorrected >= accuracyPx);
    tpr(i) = TP/(TP+FN);
    pre(i) = TP/(TP+FP);
    
    TP2 = sum(circCorrected < accuracyPx+1);
        FN2 = length(allvals) - TP2;

    tpr2(i) = TP2/(TP2+FN2);

    
end

% end
f=figure,plot(plist,tpr)
hold on
plot(plist,tpr2)
xlabel('pvalue thresh')
legend({'Exact','+-1'})
title('True positive rate')
% hold on
% plot(pre)

