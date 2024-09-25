function [w,lA]=precalc_pval_test_pcc_single_S2(output)

rng("default")

% Parameters to generate random data:
N = 20;% number of random fragments % default 200
pxpsf = 2.72; % psf
mL = 1250; % mean length of fragment % make sure long enough
stdL = 100; %std of length of fragment
snr = 2:1:10;
stdE = 0.05; %0.02; % fragment length-rescale factor std 0.02
isC = 1; % whether circular barcode
tL = 200000; % rand length in px % default 10000
minL = 150; % minimum length / same as minoverlap?
stdM = 20; % additive noise mean
sF = 0.9:0.025:1.1; % length re-scaling factor

% option: generate new for each
[bGrand, barLong, origPos, origFlip, origStr] = gen_synth_data(tL, pxpsf, ...
    N, mL, stdL, stdM, stdN, stdE, isC); % keep all

% snr that we choose here
i = 3;

stdN = sqrt(stdM.^2/snr(i)); % additive noise std. Signal mean=200 is hardcoded

import Core.gen_synth_data;


w = 100:20:400;
lA = 500:50:800;
maxPCC = cell(length(lA),length(w));
osOutput = cell(length(lA),length(w));
for j = 1:length(w)
    for k=1:length(lA)
        j
        curLen = lA(k);
        minOverlap = w(j);
        %%


        % cut bGrand to wanted length:
        barcodes = cellfun(@(x) x.rawBarcode(x.rawBitmask),bGrand,'un',false);
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
%%
save(output,'maxPCC','osOutput','lA','w');

end
