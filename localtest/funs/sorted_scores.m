
function [sortedVals, sortedIds,pscores,fullScores,overlaplen,partialScore,partialLength,pval,pvalLeftOver,pvalCombined] = sorted_scores(oS)
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

    psc = reshape([oS(:).score], size(oS,1),size(oS,2)); % from bg_test_1
    filtM = triu(psc);
    filtM(filtM==0) = nan;
    pscores = filtM(:);
    [sortedVals, sortedIds] = sort(pscores,'desc','MissingPlacement','last');

    fullScores = reshape([oS(:).fullscore], size(oS,1),size(oS,2)); % from bg_test_1
    fullScores = triu(fullScores);
    fullScores(fullScores==0) = nan;
    fullScores = fullScores(:);
    fullScores = fullScores(sortedIds);

    overlaplen = reshape([oS(:).overlaplen], size(oS,1),size(oS,2)); % from bg_test_1
    overlaplen = triu(overlaplen);
    overlaplen(overlaplen==0) = nan;
    overlaplen = overlaplen(:);
    overlaplen = overlaplen(sortedIds);

    if nargout > 5
        partialScore = reshape([oS(:).partialScore], size(oS,1),size(oS,2)); % from bg_test_1
        partialScore = triu(partialScore);
        partialScore(partialScore==0) = nan;
        partialScore = partialScore(:);
        partialScore = partialScore(sortedIds);
        %   
        partialLength = reshape([oS(:).partialLength], size(oS,1),size(oS,2)); % from bg_test_1
        partialLength = triu(partialLength);
        partialLength(partialLength==0) = nan;
        partialLength = partialLength(:);
        partialLength = partialLength(sortedIds);
    end


    % some p-value for sorted scores
    if nargout > 7
        w = oS(1,2).h;
        lenA = reshape([oS(:).lenA], size(oS,1),size(oS,2)); % from bg_test_1
        lenA = triu(lenA);
        lenA(lenA==0) = nan;
        lenA = lenA(:);
        lenA = lenA(sortedIds);

        lenB = reshape([oS(:).lenB], size(oS,1),size(oS,2)); % from bg_test_1
        lenB = triu(lenB);
        lenB(lenB==0) = nan;
        lenB = lenB(:);
        lenB = lenB(sortedIds);

        % p-value for full
        import Zeromodel.beta_ev_cdf;
        pvalfun = @(x,l1,l2,nuF) 1-beta_ev_cdf(x,nuF*w,1,2*max(l2,l1),0);

        k=length(sortedVals);
        minL = 50;
        pval = nan(1,k);

        for i=1:k
            pval(i) = pvalfun(sortedVals(i),lenA(i),lenB(i),0.07);
        end
        %p-pvalue for single
        import Zeromodel.beta_cdf;
        pvalfunLeftover = @(x,l1,nuF) 1-beta_cdf(x,nuF*l1,1,0);
        pvalLeftOver = nan(1,k);
        for i=1:k
            if partialLength(i) > minL
                pvalLeftOver(i) = pvalfunLeftover(partialScore(i),partialLength(i),0.07);
            end
        end
%         pval = [];
        pvalCombined = (norminv(1-pval)+norminv(1-pvalLeftOver))/sqrt(2);
        pvalCombined(isnan(pvalCombined)) = (norminv(1-pval(isnan(pvalCombined))));
    end
end
