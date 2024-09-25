function [ theoryCurveUnscaled_pxRes, bitmask, probSeq,theorSeq] = compute_hca_theory_barcode( seq,sets, overlapLength )
    % compute_hca_theory_barcode Computes theory barcode for hca
    %
    % :param seq: sequence (oriented at human chromosomes)
    % :param sets: settings file
    % :param overlapLength: overlap length, which should be quite long
    % :returns: theoryCurveUnscaled_pxRes, bitmask, probSeq,theorSeq
    % 
    % TODO: to be removed if compute_theory_barcode produces the same
    % results for fixed overlap/etc parameters given here
    
    import CBT.Hca.Core.Theory.cb_theory;

    if nargin < 3
        overlapLength = 10000;
    end
    
    % find bp positions where the letters are undefined or unknown, these
    % will be converted to bitmask later
    undefinedBasepairs = ones(1,length(seq));
    positions = find(seq=='N');
    undefinedBasepairs(positions) = zeros(1,length(positions));
    
    % We do not want to include regions with a lot of N's, so we randomize
    % these parts. Later we could add
    if sets.theoryGen.isLinearTF
        seq = [repmat('A',1,10000) seq repmat('A',1,10000)];
    end
    
    % The only difference between 3.8.4 and 4.0:
    % convert to int to keep the same structure as HCA 4.0
    seq = nt2int(seq, 'ACGTOnly',1);

    seq(seq==0) = 15;
    % change the bp with unknown letters into random letters (for better
    % barcodes close to unknown positions
    %seq(positions) = randseq(length(positions)); % change the unknowns into random
    
    if length(seq)< 500000
        overlapLength = 0;
        lenDiv = length(seq);
    else
        lenDiv = 200000;
    end
    % we divide sequence into parts of equal length
    seqSet = cell(1,ceil((size(seq,2)-overlapLength)/lenDiv));
    
    seqSet{1} = [seq(end-overlapLength+1:end) seq(1:lenDiv+overlapLength)];


    for i=2:(size(seq,2)-overlapLength)/lenDiv
        seqSet{i}=seq(lenDiv*(i-1)-overlapLength+1:lenDiv*(i)+overlapLength);
    end
    
    if overlapLength ~= 0
        seqSet{i+1} = [seq(lenDiv*(i)-overlapLength+1:end) seq(1:overlapLength)];
    end
    
    cN = sets.theoryGen.concN;
    cY = sets.theoryGen.concY;
    untrustedRegion = 1000;
  
    values = sets.model.netropsinBindingConstant;
    yoyo1BindingConstant = sets.model.yoyoBindingConstant;
    numSeqs = length(seqSet);
    theoryProb_bpRes = cell(numSeqs, 1);

 
    for seqNum = 1:numSeqs
        ntSeq = seqSet{seqNum};
        % compute YOYO-1 binding probabilities
        probsBinding = cb_theory(ntSeq, cN,  cY,yoyo1BindingConstant, values, untrustedRegion);

        % YOYO-1 binding probabilities
        theoryProb_bpRes{seqNum} = probsBinding;
    end

    clear seqSet
 
    psfSigmaWidth_bps = sets.theoryGen.psfSigmaWidth_nm / sets.theoryGen.meanBpExt_nm;

    import CBT.Core.convert_bpRes_to_pxRes;
    meanBpExt_pixels = sets.theoryGen.meanBpExt_nm / sets.theoryGen.pixelWidth_nm;
   
	% Convolve with a Gaussian
    ker = fftshift(images.internal.createGaussianKernel(psfSigmaWidth_bps, length(theoryProb_bpRes{1})));
    multF=conj(fft(ker));
    
    probSeq = zeros(1,length(seq));
    
    theorSeq = zeros(1,length(seq));

    for i=1:length(theoryProb_bpRes)-1
        theoryBar_bpRes = ifft(fft(theoryProb_bpRes{i}).*multF); 
        probSeq(1+(i-1)*lenDiv:(i)*lenDiv) = theoryBar_bpRes(overlapLength+1:end-overlapLength);
        theorSeq(1+(i-1)*lenDiv:(i)*lenDiv) = theoryProb_bpRes{i}(overlapLength+1:end-overlapLength);
    end
    
	% Convolve with a Gaussian
    ker = fftshift(images.internal.createGaussianKernel(psfSigmaWidth_bps, length(theoryProb_bpRes{end})));
    multF=conj(fft(ker));
    
    theoryBar_bpRes = ifft(fft(theoryProb_bpRes{end}).*multF); 
    probSeq((length(theoryProb_bpRes)-1)*lenDiv+1:end) = theoryBar_bpRes(overlapLength+1:end-overlapLength);
    theorSeq((length(theoryProb_bpRes)-1)*lenDiv+1:end) = theoryProb_bpRes{end}(overlapLength+1:end-overlapLength);

    if sets.theoryGen.isLinearTF
    	probSeq = probSeq(10001:end-10000);
        theorSeq = theorSeq(10001:end-10000);
    end
    
    % pixel resoution barcode
    theoryCurveUnscaled_pxRes = convert_bpRes_to_pxRes(probSeq, meanBpExt_pixels, sets.theoryGen.isLinearTF);
    
    % todo: make this more accurate by allowing a % (1%?) of unknown
    % letters
    bitmask = undefinedBasepairs(round(1:1/meanBpExt_pixels:length(undefinedBasepairs)));
end
