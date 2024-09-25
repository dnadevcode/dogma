function [barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen, synth,PSF_WIDTH_PIXELS,NUM_RAND_FRAGMENTS)

    
% idx = 4;
% numF = 10;

if nargin < 6
    NUM_RAND_FRAGMENTS = 100;
end

if synth
%     NUM_RAND_FRAGMENTS = 100;
    if nargin <5
        PSF_WIDTH_PIXELS = 300/440;
    end
%     RAND_LENGTH_MIN = 1500;
    % RAND_LENGTH_2 = 800;
    
    % sF = 1;%
    
    % sF = 0.8:0.01:1.2;
    
    import Nullmodel.gen_random;
    [~,barcodeGen1] = gen_random(NUM_RAND_FRAGMENTS/2,PSF_WIDTH_PIXELS,numF,0);
    [~,barcodeGen2] = gen_random(NUM_RAND_FRAGMENTS/2,PSF_WIDTH_PIXELS,minLen,0);

%
else
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    
    
    import Core.create_barcodegen;
    [barcodeGen1, kymoStructs1, sets1] = create_barcodegen(testSet{1},1,numF,1,timestamp,1);
    if length(testSet)==2
        [barcodeGen2, kymoStructs2, sets2] = create_barcodegen(testSet{2},1,numF,1,timestamp,1);
    end
end
%%
% only keep barcodes longer than minimum length.
% minLen = 300; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen1);
barcodeGen1 = barcodeGen1(barLens>=minLen);
lengths1 = cellfun(@(x) sum(x.rawBitmask),barcodeGen1);

if length(testSet)==2 || synth
    barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen2);
    barcodeGen2 = barcodeGen2(barLens>=minLen);
    lengths2 = cellfun(@(x) sum(x.rawBitmask),barcodeGen2);
else
    barcodeGen2 = [];
    lengths2 = [];
end

