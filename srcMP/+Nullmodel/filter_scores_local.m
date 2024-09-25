function [passLocalThresh,coefDiffsMP] = filter_scores_local(sortedVals,mpMaxLenBased,minOverlap,minLen,N)
    %   filter_scores_local - calculate which scores pass the local score
    %   threshold
    %
    %   Args:
    %       sortedVals - sorted local scores
    %       mpMaxLenBased - list of local mp scores
    %       minOverlap - overlap length for local scores
    %       minLen - vector of lengths for mpMaxLenBased
    %   Returns:
    %       passLocalThresh - binary vector with 1 if local length thresh
    %       passed

    if nargin < 5
        N = length(sortedVals);
    end

    N = min(N,sum(~isnan(sortedVals)));

    idxThreshMP = arrayfun(@(x) find(x<=minLen,1,'first'),minOverlap);
    coefDiffsMP = sortedVals(1:N)-mpMaxLenBased(idxThreshMP); 
    passLocalThresh = coefDiffsMP > 0;


end

