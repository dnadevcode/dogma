function [passGlobalThresh] = filter_scores_global(foS,mpMaxLenBased,overlaplen,minLen,N)
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
        N = sum(~isnan(foS));
    end

    N = min(N,sum(~isnan(foS)));


    idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),overlaplen(1:N));

    passGlobalThresh = foS(1:N) > mpMaxLenBased(idxThresh)';

end

