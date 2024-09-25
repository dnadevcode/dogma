function [sortedVals,sortedIds,fullScores,posShiftSelected] = filter_good_scores(oS,mpMax,minLen,thresCC,posShift)
%   filter_good_scores
%
%   Args:
%       oS - pair-overlap result
%       mpMax - mp score thresholds
%       minLen - minimum length
%       thresCC - threshold for local scores
%       posShift - shifts from the GT position

%   Returns:
%       sortedVals - sorted local scores
%       sortedIds - pair indexes for sorted scores
%       fullScores - full scores sorted by local scores
%       lastCC - number of local scores below thresh
%       
[sortedVals, sortedIds, pscores, fullScores, overlaplen] = sorted_scores(oS);


[lastCC,~] = find(sortedVals<thresCC,1,'first');
fullScores = fullScores(1:lastCC);
sortedVals = sortedVals(1:lastCC);

%thresh for each length
idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),overlaplen(1:lastCC));

idxs = fullScores>mpMax(idxThresh)';
sortedVals = sortedVals(idxs);
sortedIds = sortedIds(idxs);

fullScores = fullScores(idxs);

if nargin > 4
    posShiftSelected = posShift(idxs);
else
    posShiftSelected = [];
end

end

