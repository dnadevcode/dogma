function [data] = compute_random_max_cc(len1, rand2, psf, numRnd)
    % compute_random_max_cc

    % Computes maximum correlation for a given number of random barcodes

    % Pregenerates p-value database file
    %     Args:
    %         len1,rand2,psf, numRnd
    % 
    %     Returns:
    %         data  
    %     Example:
    
    bit1 = ones(1,len1);
    
    data = zeros(1,numRnd);
    import CBT.Hca.Core.Pvalue.convolve_bar;
	import SignalRegistration.unmasked_pcc_corr;

    for i =1:numRnd
        rand1 = normrnd(0,1, 1, len1);
        rand1 = convolve_bar(rand1, psf, len1 );
        [ccM] = unmasked_pcc_corr(rand1, rand2, bit1);
        data(i) = max(ccM(:));
    end
    
end

