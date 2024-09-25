function [ccthresh] = ccthresh_for_comparisonstruct(minLen,nuF,lengthsThry, pvalThresh)
    % cc-thresh based on different barcode lengths. Global score
    %
    %   Args:
    %       minLen - lengths
    %       nuF - evd parameter (depends on psf)
    %       lengthsThry - approximate length theory
    %
    %   Returns:
    %       ccthresh

    if nargin < 4
        pvalThresh = 0.01;
    end

    import Zeromodel.beta_ev_cdf;

    % form for Neff (instead of fitting EVD). Check pval_test_loop tests for argumentation of
    % correct form / non-circ
    funNeffective = @(lengthsThry) 2*(lengthsThry);

    ccthresh = zeros(1,length(minLen));
    % loop
    for i=1:length(minLen)
        pvalfun = @(x) 1-beta_ev_cdf(x,nuF*minLen(i),1,funNeffective(lengthsThry),0)-pvalThresh;
        ccthresh(i) = fzero(pvalfun,0.5);
    end

end

