function [pccs,lengths,parameters] = gen_scores_test() % choose some params

    % first generate 1000 random scores
    N = 100;% number of random fragments
    pxpsf = 2.72; % psf
    mL = 850; % mean length of fragment
    stdL = 0; %std of length of fragment
    snr = 10;
    stdE = 0.05;%0.02; % fragment length-rescale factor std 0.02
    isC = 1; % whether circular barcode
    tL = 10000; % rand length in px
    minL = 150; % minimum length / same as minoverlap?
    stdM = 20; % additive noise mean

    sF = 0.9:0.01:1.1;

    

    stdN = sqrt(stdM.^2/snr); % additive noise std. Signal mean=200 is hardcoded

    % synthetic (not based on theory)
    import Core.gen_synth_data;
    [bG, barLongTrue, origPos, origFlip, origStr] = gen_synth_data(tL, pxpsf, ...
        N, mL, stdL, stdM, stdN, stdE, isC);

    tL2 = 1000;
    % single long
    [~, barLong, ~, ~, ~] = gen_synth_data(tL2, pxpsf, ...
        1, mL, stdL, stdM, stdN, stdE, isC);


% maxPCC=cell(1,NN);

    % values..
import mp.mp_profile_stomp_nan_dna;
comparisonFun2 = @(X1,X2, bitX1, bitX2, w, kk) mp_profile_stomp_nan_dna(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))),50);



comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),0,50);

w = 300;
pccs = zeros(1,length(bG));
for i=1:length(bG)
    i
    A = arrayfun(@(y) imresize(bG{i}.rawBarcode(bG{i}.rawBitmask),'Scale' ,[1 y]),sF,'un',false);
    pcCur = 0;
    for j=1:length(A)
%         pcCur = max([pcCur,comparisonFun(A{j},barLong,ones(1,length(A{j})),ones(1,length(barLong)),2^16)]);
        pcCur = max([pcCur,comparisonFun2(A{j},barLong,ones(1,length(A{j})),ones(1,length(barLong)),w,2^16)]);
    end
    pccs(i) = pcCur;
end


lengths = [mL,tL2,w]
import CBT.Hca.Core.Pvalue.compute_evd_params;
[ parameters ] = compute_evd_params( pccs,mL/4 );

% 
% 
% 
% for k=1:NN
%     % run loop
% %     tic
% %     k
%     maxPCC{k} = nan(1,NUM_RAND_FRAGMENTS);
%     
% 
%     [barcodes] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN+(k-1)*10,0);
%     %
%     bitmasks = cellfun(@(x) true(size(x)), [barcodes], 'un', 0);
% 
% 
%     % create bar structure 
%     barSynth = cell2struct([barcodes; bitmasks]',{'rawBarcode','rawBitmask'},2);
% 
% 
%     [barcodes2] = gen_random(1,PSF_WIDTH_PIXELS,RAND_LENGTH_2,0);
%     bitmasks2 = cellfun(@(x) true(size(x)), [barcodes2], 'un', 0);
% 
%     barSynth2 = cell2struct([barcodes2; bitmasks2]',{'rawBarcode','rawBitmask'},2);
% 
% 
%     % barSynth = [cell2struct(barcodes,{'rawBarcode'},1);cell2struct(bitmasks,{'bitmasks'},1)];
% 
%     % have to save separately..
% 
%     maxPCC{k}= pccs{k};
% end


end

