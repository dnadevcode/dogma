function [] = precalc_test_pcc_single_real_S3_S4()

% load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/2all_2023-11-28_09_43_54.mat')
load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/3all_2023-11-24_12_38_55.mat')

% S3

%%
minL = 500;

bGrand = barcodeGen;


lengths = cellfun(@(x) sum(x.rawBitmask),bGrand);

bGrand = bGrand(lengths>minL);

lengths = lengths(lengths>minL);
[maxL,idx] = max(lengths);

% as bGrand2, the longest barcode
bGrand2 = bGrand(idx);
bGrand(idx) = [];

bGrand = bGrand(500:600);
%%


rng("default")

% Parameters to generate random data:
N = 100;% number of random fragments % default 200
pxpsf = 2.72; % psf
mL = 2850; % mean length of fragment % make sure long enough
stdL = 100; %std of length of fragment
snr = 0.5:0.5:5;
stdE = 0.05; %0.02; % fragment length-rescale factor std 0.02
isC = 1; % whether circular barcode
tL = 200000; % rand length in px % default 10000
minL = 150; % minimum length / same as minoverlap?
stdM = 20; % additive noise mean
sF = 0.9:0.025:1.1; % length re-scaling factor



% w = 100:20:400;
lA = 200:20:350;
lB = 500:50:maxL;

comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),0,50);

maxPCC = cell(length(lA),length(lB));
osOutput = cell(length(lA),length(lB));
pccs =  cell(length(lA),length(lB)); 
for j = 1:length(lB)
    for k=1:length(lA)
    j
    curLen = lA(k);
    curLenB = lB(j);
    %%

    % cut bGrand to wanted length:
    barcodes = cellfun(@(x) x.rawBarcode(x.rawBitmask),bGrand','un',false);
    barcodes = cellfun(@(x) x(randi(length(x)-curLen)+(1:curLen)),barcodes,'un',false); % random cut-out of length curLen
    bitmasks = cellfun(@(x) true(size(x)), barcodes, 'un', 0);
    barSynth = cell2struct([barcodes, bitmasks],{'rawBarcode','rawBitmask'},2);

        % cut bGrand to wanted length:
    barcodes2 = cellfun(@(x) x.rawBarcode(x.rawBitmask),bGrand2,'un',false);
    barcodes2 = cellfun(@(x) x(randi(length(x)-curLenB)+(1:curLenB)),barcodes2,'un',false); % random cut-out of length curLen
    bitmasks2 = cellfun(@(x) true(size(x)), barcodes2, 'un', 0);
    barSynth2 = cell2struct([barcodes2, bitmasks2],{'rawBarcode','rawBitmask'},2);


    for l=1:length(barSynth)
        A = arrayfun(@(y) imresize(barSynth(l).rawBarcode(barSynth(l).rawBitmask),'Scale' ,[1 y]),sF,'un',false);
        pcCur = 0;
        for d=1:length(A)
                pccThis = comparisonFun(A{d},barcodes2{1},ones(1,length(A{d})),ones(1,length(barcodes2{1})),2^16);
                pcCur = max([pcCur,pccThis]);
        end
        pccs{j,k}(l) = pcCur;
    end
    end
end

save('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S1_data_updatedReal.mat','pccs','lA','lB');

%% S4

bGrandN = bGrand(1:20);
w = 100:20:300;
lA = 340:20:500;
maxPCC = cell(length(lA),length(w));
osOutput = cell(length(lA),length(w));
for j = 1:length(w)
    for k=1:length(lA)
    j
    curLen = lA(k);
    minOverlap = w(j);
    %%


    % cut bGrand to wanted length:
    barcodes = cellfun(@(x) x.rawBarcode(x.rawBitmask),bGrandN','un',false);
    barcodes = cellfun(@(x) x(randi(length(x)-curLen)+(1:curLen)),barcodes,'un',false); % random cut-out of length curLen
    bitmasks = cellfun(@(x) true(size(x)), barcodes, 'un', 0);
    barSynth = cell2struct([barcodes, bitmasks],{'rawBarcode','rawBitmask'},2);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(barSynth,sF, minOverlap, [],29);
    toc
    localScore = [oS(:).score]; % local
    lenA = [oS(:).lenA]; % lenA
    lenB = [oS(:).lenB]; % lenB\
    partialScore = [oS(:).partialScore];
    partialLength = [oS(:).partialLength];
    overlapLength = [oS(:).overlaplen];

    idx = reshape(1:size(oS,1)*size(oS,2), size(oS,1),size(oS,2)); % from bg_test_1
    idx = tril(idx);

    localScore(idx(idx~=0)) = nan;
    partialScore(idx(idx~=0)) = nan;

    maxPCC{k,j} = localScore(~isnan(localScore));
    osOutput{k,j} = oS;
    end
end

save('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S2_data_updated2Real.mat','maxPCC','osOutput','lA','w');




end

