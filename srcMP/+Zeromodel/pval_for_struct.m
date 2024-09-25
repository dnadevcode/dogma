function [pscores] = pval_for_struct(overlapStruct, ccthresh, nuF)


% % pvals
% ccthresh = 0.7;
% nuF = 0.05; % dep on PSF.


totL = length([overlapStruct.sc]);
pscores = nan(1,totL);
import Zeromodel.beta_ev_cdf

for ii=1:totL
%     ii
    if overlapStruct(ii).sc > ccthresh
        % todo: simplify for quicker calculation
        pscores(ii) =( 1-beta_ev_cdf(overlapStruct(ii).sc, nuF*overlapStruct(ii).overlaplen, 1, 2*max(overlapStruct(ii).lenA,overlapStruct(ii).lenB)));
    end
end
end

