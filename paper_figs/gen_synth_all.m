
rng('default');


N = 200;% number of random fragments
pxpsf = 2.72; % psf
mL = 850; % mean length of fragment
stdL = 100; %std of length of fragment
snr = 1.5;
stdE = 0.05;%0.02; % fragment length-rescale factor std 0.02
isC = 1; % whether circular barcode
tL = 10000; % rand length in px
minL = 150; % minimum length / same as minoverlap?

% synthetic data (if not generated yet
import Core.load_synth_data;
[bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(N,2.72,850,100,snr,0.05,1,10000,150);

mean(cellfun(@(x) x.maxcoef,synthStr{1}))

% 
% % synthetic data (mostly for evaluation)
% import Core.load_synth_data;
% [bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data();

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sF = 0.95:0.01:1.05; % length re-scaling factor
minOverlap = 300; % minimum overlap

resRun = cell(1,length(bG));
for idxRun = 1:length(bG)
    barcodeGen = bG{idxRun};
    lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(barcodeGen',sF, minOverlap, timestamp);
    toc
    resRun{idxRun}.oS = oS;
end
save(['synth','_',timestamp,'resRun.mat'] ,'resRun','bG','synthStr2','synthStr','theoryStruct','mpMaxLenBased','minLen' )

