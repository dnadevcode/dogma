function [ theoryCurveUnscaled_pxRes, bitmask ] = compute_theory_barcode( seq,sets )
    % Computes theory barcode for hca
    
    % input 
    % seq - sequence (oriented at human chromosomes)
    % sets - settings
    
    % output 
    % barcode
    
    % We do not want to include regions with a lot of N's, so we randomize
    % these parts. Later we could add
%     if sets.isLinearTF
%         seq = [randseq(1000) seq randseq(1000)];
%     end
    
    seq( find(seq=='N')) = randseq(length(find(seq=='N'))); % change the unknowns into random
    
    concNetropsin_molar = sets.concNetropsin_molar;
    concYOYO1_molar = sets.concYOYO1_molar;

    import CBT.Core.cb_netropsin_vs_yoyo1_plasmid;
    probsBinding = cb_netropsin_vs_yoyo1_plasmid(seq, concNetropsin_molar,  concYOYO1_molar);

    % YOYO-1 binding probabilities
    theoryProb_bpRes = probsBinding.Yoyo1;
        

    import Microscopy.Simulate.Core.apply_point_spread_function;
    psfSigmaWidth_bps = sets.psfSigmaWidth_nm / sets.meanBpExt_nm;

    import CBT.Core.convert_bpRes_to_pxRes;
    meanBpExt_pixels = sets.meanBpExt_nm / sets.pixelWidth_nm;

%     isLinearTF = sets.isLinearTF;
%     widthSigmasFromMean = sets.widthSigmasFromMean;

%  psfSigmaWidth_bps = 300;
    import CA.CombAuc.Core.Zeromodel.gaussian_kernel;
    ker = gaussian_kernel(length(theoryProb_bpRes), psfSigmaWidth_bps);
    multF=conj(fft(ker));

    probSeq = ifft(fft(theoryProb_bpRes).*multF); 
 
%     if sets.isLinearTF
%     	probSeq = probSeq(1001:end-1000);
%     end
    
    theoryCurveUnscaled_pxRes = convert_bpRes_to_pxRes(probSeq, meanBpExt_pixels);

    bitmask = ones(1,length(theoryCurveUnscaled_pxRes));
end

