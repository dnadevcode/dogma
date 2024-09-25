figs% validate ground-truth

% Generate simulated barcodes for given SNR
% Compare to theory
% Calculate MP length based statistics
% Report how many end up at correct place


% Over different signal-to-noise-ratio
% addpath(genpath('/home/etlar/albertas/reps/hca'))
% addpath(genpath('/home/etlar/albertas/reps/bargroupingprototype'))
% addpath(genpath('/home/etlar/albertas/reps/discriminative_sequences'))


timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%%


allowedDist = 20;
NUM_RAND_FRAGMENTS = 50;% number of random fragments
PSF_WIDTH_PIXELS = 2.72; % psf
MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
FRAGMENT_LENGTH_STD = 100; %std of length of fragment
ADDED_NOISE_MEAN = 0.4; % additive noise mean

FRAGMENT_STRETCH_STD = 0.05;%0.02; % fragment length-rescale factor std 0.02
IS_CIRCULAR = 1; % whether circular barcode
TOTAL_RAND_LENGTH = 10000; % rand length in px
minL = 100; % minimum length / same as minoverlap?


%% min length calculation for different overlap lengths
nmPerPx = 110;
nmbp = 0.2;
minLen =[minL:50:3000]; % min to max
[MP,mpMaxLenBased,theoryStructRev,MPI,~] = bargrouping_minimum_length([],nmPerPx,nmbp,1,30,minLen, TOTAL_RAND_LENGTH*nmPerPx/nmbp);
minOverlapList = 100:50:800; % minimum overlap
sF = 0.85:0.025:1.15; % length re-scaling factor

%%
%
SNR = 8;
tpr = zeros(1,length(SNR));
fpr = zeros(1,length(SNR));

cdiff = 0.05;

for i=1:length(SNR)

    ADDED_NOISE_STD = 10/SNR(i); % additive noise std. Signal mean=10 is hardcoded
    
    %% generate synthetic data
    import Nullmodel.gen_rand;
    [barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr] = gen_rand(...
        NUM_RAND_FRAGMENTS, PSF_WIDTH_PIXELS,MEAN_FRAGMENT_LENGTH,FRAGMENT_LENGTH_STD,...
        ADDED_NOISE_MEAN,ADDED_NOISE_STD,FRAGMENT_STRETCH_STD,IS_CIRCULAR,minL,TOTAL_RAND_LENGTH);%,barcodeGenT{1}.rawBarcode);

    barcodeGen = barcodeGen';

    goodBars = zeros(1,length(minOverlapList));
    for idxRun = 1:length(minOverlapList)
        minOverlap = minOverlapList(idxRun);
    
        lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
        tic
        % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
        [oS] = calc_overlap_mp(barcodeGen(lengths>=minOverlap),sF, minOverlap,timestamp);
        toc
        import Core.plot_match_pcc;
        [sortedVals, sortedIdsAll,pscores] = sorted_scores(oS);
        idxThreshMP = arrayfun(@(x) find(x<=minLen,1,'first'),minOverlap);
        coefDiffsMP = sortedVals-mpMaxLenBased(idxThreshMP);
        barsPassThreshMP = find(coefDiffsMP > 0);
        
        thresCC = mpMaxLenBased(idxThreshMP);
        import Core.filter_good_scores
        [sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBased,minLen,thresCC);
        
        goodBars(idxRun) = length(sortedVals);
    end

end

f = figure
plot(minOverlapList,goodBars)
hold on
plot([300 300],[0 max(goodBars)],'red-')
xlabel('Min overlap length');
ylabel('Number of significant overlaps')
legend({'Number of significant overlaps','overlap length w'},'location','southoutside')
    print('FIGS/FigS11.eps','-depsc','-r300');

%     
% 
%     sets.comparisonMethod = 'mass_pcc';
%     
%     sF = 0.9:0.01:1.1;
%     [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(barcodeGen,theoryStruct,sF,sets);
%     
%     allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
%     allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);
%     
%     idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
%     coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
%     
%     barsParsThresh = coefDiffs > cdiff;
    % barsPassThresh = find();

    
    % super_quick_plot(6,barcodeGen,comparisonStruct,theoryStruct)
    
    %% check success rate statistics
    N = length(barcodeGen);
    posShift = zeros(1,N);
    tp = zeros(1,N);
    fn = zeros(1,N);
    fp = zeros(1,N);
    tn = zeros(1,N);
    
    success = 0;
    for idx=1:N
        pos = comparisonStruct{idx}.pos(1);
        posGT = synthStr{idx}.pos;
        %         posDifTheory = pos-posGT;
        posShift(idx) = min([dist(pos,posGT) dist(pos+TOTAL_RAND_LENGTH,posGT) dist(pos-TOTAL_RAND_LENGTH,posGT) ]);
        if posShift(idx) <= allowedDist
            if barsParsThresh(idx)==1
                tp(idx) = 1;
            else
                fn(idx) = 1;
            end
        else
            if barsParsThresh(idx)==1
                fp(idx) = 1;
            else
             tn(idx) = 1;
            end
        end
    end
    
    tpr(i) = sum(tp)/length(barcodeGen);%(sum(tp)+sum(fn)); % sensitivity
    fpr(i) = sum(fp)/(sum(fp)+sum(tn)); % fpr
end

figure,plot(tpr)
hold on
plot(fpr)
xlabel('SNR')
ylabel('Success rate')
legend({'tpr','fpr'},'location','southoutside')
title('Validation of ground truth method')
print('FIGS/FigS10.eps','-depsc','-r300');


%% Validation 2