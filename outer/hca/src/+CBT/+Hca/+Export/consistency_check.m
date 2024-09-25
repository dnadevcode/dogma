function [] = consistency_check( b1,cutB,dd,bitm)
    if nargin < 4
        bitm = ones(1,length(b1));
    end
 
    % Consistency check that the cc for the cut out barcodes are the same
    import SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs;
    [xcorrs, ~, ~] = get_no_crop_lin_circ_xcorrs(b1,cutB', bitm, ones(1,length(cutB)));

    % d1 = Comparison.cc_fft(zscore(b1),zscore(cutB'));
    
    xcorrs(1,1)
    
    if (abs(xcorrs(1,1)-dd) < 10^-10)
        display('Theory vs exp. plot passes consistency check');
    else
        display('Theory vs exp. plot is not consistent');
    end

end

