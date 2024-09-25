function []=precalc_pval_test_pcc_single_simple_S1(output)
%% Generate synthetic data
rng("default")

% Parameters to generate random data:
N = 100;% number of random fragments % default 200
pxpsf = 2.72; % psf
mL = 2850; % mean length of fragment % make sure long enough
stdL = 100; %std of length of fragment
snr = 2:1:10;
stdE = 0.05; %0.02; % fragment length-rescale factor std 0.02
isC = 1; % whether circular barcode
tL = 200000; % rand length in px % default 10000
minL = 150; % minimum length / same as minoverlap?
stdM = 20; % additive noise mean
sF = 0.9:0.025:1.1; % length re-scaling factor

i = 3;
stdN = sqrt(stdM.^2/snr(i)); % additive noise std. Signal mean=200 is hardcoded

import Core.gen_synth_data;

% option: generate new for each
[bGrand, barLong, origPos, origFlip, origStr] = gen_synth_data(tL, pxpsf, ...
    N, mL, stdL, stdM, stdN, stdE, isC); % keep all
[bGrand2, barLong, origPos, origFlip, origStr] = gen_synth_data(tL, pxpsf, ...
    1, mL, stdL, stdM, stdN, stdE, isC); % keep all


%% Run comparison

% w = 100:20:400;
lA = 300:20:700;
lB = 800:50:2000;

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
        barcodes = cellfun(@(x) x.rawBarcode(x.rawBitmask),bGrand,'un',false);
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

%% save comparison to output
save(output,'pccs','lA','lB');

end
