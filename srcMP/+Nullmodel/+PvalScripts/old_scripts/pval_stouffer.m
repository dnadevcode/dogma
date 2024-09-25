% create Stouffer histogram for p-values local/partial

import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 100;
% minOverlapPx = 300;
PSF_WIDTH_PIXELS = 300/110;
randLenPx = 5000;
randLenPx2 = 5000;

[~ , barSynth] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,randLenPx2,0);
% [~, barSynth2] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,randLenPx,0);

sF = 0.9:0.01:1.1;
minOverlap = 300;
timestamp = 'temp';

[oS] = calc_overlap_mp(barSynth,sF, minOverlap, timestamp);


% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver,pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,minOverlap, 0.08,inf); %localCCStruct

