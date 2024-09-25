
function [sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS, w, alphaNu1,pthresh,localCCStruct,alphaN1,alphaNu2,alphaN2)
    % sorted_scores - sorted alignment scores  
    %   Args:
    %       oS - overlap result structure
    %
    %   
    %   Returns:
    %       sortedVals - sorted local scores
    %       sortedIds - sorted ids for local scores
    %       pscores - ? 
    %       fullScores - full overlap scores
    %       overlaplen - overlap length
    %       partialScore - partial (leftover) score
    %       partialLength - partial (leftover) length
    %       pval - pvalue local
    %       pvalLeftOver - pvalue leftover
    %       pvalCombined - pvalue combined (Stouffer's method)
    %       sortedValsBad - sorde

    % todo: should also work for global-overlaps, i.e. w = nan;

    % get all the scores
    localScore = [oS(:).score]; % local
    lenA = [oS(:).lenA]; % lenA
    lenB = [oS(:).lenB]; % lenB\
    partialScore = [oS(:).partialScore];
    partialLength = [oS(:).partialLength];
    overlapLength = [oS(:).overlaplen];

    idx = reshape(1:size(oS,1)*size(oS,2), size(oS,1),size(oS,2)); % from bg_test_1
    idx = tril(idx);

    localScore(idx(idx~=0)) = nan;
    partialScore(idx(idx~=0)) = nan;

    % Alternative method using pre-generated list of CC thresholds for
    % specific w
    if nargin >=5 && ~isempty(localCCStruct)
        idxThreshMP = arrayfun(@(x) find(x<=localCCStruct.minLen,1,'first'),w);
        coefDiffsMP = localScore-localCCStruct.maxBasedCC(idxThreshMP); % has to be already calculated
        barsFailCCThresh = find(coefDiffsMP <= 0);
        localScore(barsFailCCThresh) = nan;
    end

    % pval MP
    if nargin < 6
        alphaN2 = 0.001; %should be rather insensitive to specific value, check figure S5-S6;
        alphaN1 =  0.2; %should be rather insensitive to specific value, check figure S5-S6;
        alphaNu2 = 0.14;
    end

    import Zeromodel.beta_ev_cdf; % correct form?
%     pvalfun = @(x,l1,l2,nuF,w) 1-beta_ev_cdf(x,nuF*w,1,nuF*2*(max(l1,l2)-w+1),0);
    pvalfun = @(x,l1,l2,nuF,w) 1-beta_ev_cdf(x,nuF*w,1,max(alphaN1*2*(max(l1,l2)-w),alphaN2*2*(l1-w+1)*(l2-w+1)),0);


    % p-pvalue for leftover
    import Zeromodel.beta_cdf;
    pvalfunLeftover = @(x,l1,nuF) 1-beta_cdf(x,nuF*l1,1,0);


    k = length(localScore);
    minL = 50; % don't calculate left-over score when length less than minL

    % p-values
    pvalLocal = nan(1,k);
    pvalLeftOver = nan(1,k);

    for i=1:k
        if isnan(w) 
            if ~isnan(overlapLength(i)) && overlapLength(i) > minL
                pvalLocal(i) = pvalfun(localScore(i),lenA(i),lenB(i),alphaNu2,overlapLength(i));
            else
                pvalLocal(i) = nan;
            end

        else
            pvalLocal(i) = pvalfun(localScore(i),lenA(i),lenB(i),alphaNu2,w);
        end

        if partialLength(i) > minL
%             pvalLeftOver(i) = pvalfun(localScore(i),lenA(i),lenB(i),nuF);
            pvalLeftOver(i) = pvalfunLeftover(partialScore(i),partialLength(i),alphaNu1);
        end
    end

    pvalCombined = (norminv(1-pvalLocal)+norminv(1-pvalLeftOver))/sqrt(2);
    pvalCombined(isnan(pvalCombined)) = (norminv(1-pvalLocal(isnan(pvalCombined))));

    pvalAll = pvalCombined;

    idxsAboveThreshLeftover = logical(~isnan(pvalLeftOver).*[pvalLeftOver>pthresh]);
    idxsAboveThreshLocal = logical(~isnan(pvalLocal).*[pvalLocal>pthresh]);
    % this completely removes values not passing threshold
    pvalCombined(idxsAboveThreshLeftover) = nan; % remove values not passing the threshold
    pvalCombined(idxsAboveThreshLocal) = nan; % remove values not passing the threshold

    pvalAll(and(~idxsAboveThreshLeftover,~idxsAboveThreshLocal)) = nan;
%     pvalAll(~idxsAboveThreshLeftover) = nan;


    [sortedVals, sortedIds] = sort(pvalCombined,'desc','MissingPlacement','last');

%     sortedVals(sortedIds(isnan(pvalCombined))) = nan;

    localScore = localScore(sortedIds);
    partialScore = partialScore(sortedIds);
    lenA = lenA(sortedIds);
    lenB = lenB(sortedIds);
    partialLength = partialLength(sortedIds);
    pvalLocal = pvalLocal(sortedIds);
    pvalLeftOver = pvalLeftOver(sortedIds);

%     idxBad = 
    [sortedValsBad, sortedIdsBad] = sort(pvalAll,'desc','MissingPlacement','last');


    sortedIdsBad = sortedIdsBad(~isnan(sortedValsBad));
    sortedValsBad = sortedValsBad(~isnan(sortedValsBad));


    %% Local scores
%     stoufMP = norminv(1-pvalLocal);
%     stoufLeft = norminv(1-pvalLeftOver);
% 
%     figure,histogram(stoufMP);hold on
%     histogram(stoufLeft)

end
