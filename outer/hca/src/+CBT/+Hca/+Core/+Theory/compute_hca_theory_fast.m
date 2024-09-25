function [ theoryCurveUnscaled_pxRes] = compute_hca_theory_fast( chr1,sets)
    % compute_hca_theory_barcode Computes theory barcode for hca
    %
    %   Args:
    %        chr1,sets, overlapLength
    %   Returns:
    %        theoryCurveUnscaled_pxRes, bitmask, probSeq, theorSeq
    % :param seq: sequence (oriented at human chromosomes)
    % :param sets: settings file
    % :param overlapLength: overlap length, which should be quite long
    % :returns: theoryCurveUnscaled_pxRes, bitmask, probSeq,theorSeq
    
	% use reproducible random numbers. Fix them for the whole calculation,
	% so that we wouldn't have strange periodic structures in the barcode
    %

    
    % overlap length for the calculation of theory
    sets.overlapLength = 3000;
    sets.blockSize = 500000;
    sets.untrustedRegion = 1000;
    
    

        
    % total length
    totalLen = numel(chr1.Data);
    
    % divide data into a number of blocks
    numNT = numel(chr1.Data);
    blockSize =  sets.blockSize;
    overlapLength =  sets.overlapLength;
    numBlocks = floor(numNT/blockSize);
    
%     % free concentrations of yoyo and netropsin
%     cN = sets.theoryGen.concN;
%     cY = sets.theoryGen.concY;
%     
%     % untrusted region
%     untrustedRegion =  sets.untrustedRegion;
%   
%     % binding probabilities for the competitive binding model
%     values = sets.model.netropsinBindingConstant;
%     yoyo1BindingConstant = sets.model.yoyoBindingConstant;

    % psf in basepairs, we need to convolve with such Gaussian
    psfSigmaWidth_bps = sets.theoryGen.psfSigmaWidth_nm / sets.theoryGen.meanBpExt_nm;

    % This function converts bpRes to pxRes
    import CBT.Core.convert_bpRes_to_pxRes;
    meanBpExt_pixels = sets.theoryGen.meanBpExt_nm / sets.theoryGen.pixelWidth_nm;

%     numSeqs = length(numBlocks+1);
%     theoryProb_bpRes = cell(numSeqs, 1);

    disp('started generating theory barcode.. 0%')
    tic
    probSeq = zeros(1,totalLen);
    
    for count=1:numBlocks+1
        % calculate the indices for the block
        start = 1 + blockSize*(count-1);
        stop = blockSize*count;
        
        % extract the block. Because of the first and last vectors becomes
        % a little bit more difficult
        if count == 1
            if sets.theoryGen.isLinearTF
                block = [ones(1,10000) ones(1,overlapLength) chr1.Data(start:stop+overlapLength)'];
            else
                block = [chr1.Data(end-overlapLength+1:end)' chr1.Data(1:stop+overlapLength)'];
            end
        else
            if count == numBlocks+1
                if sets.theoryGen.isLinearTF
                    block = [chr1.Data(start-overlapLength:end)' ones(1,overlapLength)  ones(1,10000) ];
                else
                    block = [chr1.Data(start-overlapLength:end)' chr1.Data(1:overlapLength)'];
                end
            else
                block = chr1.Data(start-overlapLength:stop+overlapLength)';
            end
        end
        
        s = rng;
        rng(0,'twister');
        block(block  > 4) = randi(4,1,sum(block > 4));
        rng(s);

    
        numWsCumSum = cumsum((block == 1)  | (block == 4) );
        probsBinding = [1 * (numWsCumSum(5:end) == numWsCumSum(1:end-4) + 4)'; 0 ;0; 0; 0];

        % Convolve with a Gaussian. The length of this differs, but only
        % three cases, which we can cover before the loop.
        ker = fftshift(images.internal.createGaussianKernel(psfSigmaWidth_bps, length(probsBinding)));
        multF = conj(fft(ker));
        
        % current theoryBar_bpRes
        theoryBar_bpRes = ifft(fft(probsBinding).*multF); 
        % now we could convert this to px resolution, but we can't because
        % overlap won't be at the exactly same places. TODO: make the
        % overlaps of exact size (so instead of 3000 bp something more accurate)
        % , so that the overlaps would match
        
        if count==numBlocks+1
            probSeq(1+(count-1)*blockSize:end) = theoryBar_bpRes(overlapLength+1:end-overlapLength);
        else
            probSeq(1+(count-1)*blockSize:(count)*blockSize) = theoryBar_bpRes(overlapLength+1:end-overlapLength);
        end
    end
    
    disp('done generating theory barcode... 100%')
    toc   
    
    if sets.theoryGen.isLinearTF
        probSeq = probSeq(10001:end-10000);
    end
        
   
    % pixel resoution barcode
    theoryCurveUnscaled_pxRes = convert_bpRes_to_pxRes(probSeq, meanBpExt_pixels, sets.theoryGen.isLinearTF);
    
    clear probSeq;
    
    % todo: make this more accurate by allowing a % (1%?) of unknown
    % letters
    %bitmask = undefinedBasepairs(round(1:1/meanBpExt_pixels:length(undefinedBasepairs)));
    

end

