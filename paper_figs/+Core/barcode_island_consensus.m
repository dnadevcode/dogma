function [aaBarcodeIslandBlockRep, Z , allscores,clusterBarcodes, rawBarcodeIslandBlockRep] = barcode_island_consensus(barcodeGen,cGenAll, ix, wmin, thresh,alphaNu)

% barcode_island_consensus - consensus with correct barcode scaling along
% overlapping windows. Uses hierarchical clustering approach. See
% consensus_real_data.mlx for example of how to use it
%

%   Args:
%       barcodeGen - all barcodes
%       cGenAll  - consensus output struct - processed table representation
%       ix - which island to calculate
%       wmin - minimum size for overlap window

%   Returns:
%       amplitudeAdjustedConsensus - amplitude adjusted consensus
%       Z - linkage
%       allscores - pairwise scores for all barcodes
%       clusterBarcodes - barcode idxs that pass the cluster threshold
%       rawConsensus - raw consensus, when no amplitude adjustment is done

% This version does single linkage algorithm. 

if nargin < 4
    wmin = 300;
end

if nargin < 5 || isempty(thresh)
    % only keep barcodes passing thresh similarity level
    thresh = 0.5;%norminv(1-0.05);
end

if nargin < 6
    alphaNu = 0.09;
end

% cgenid = cGenAll{ix}.idx;
bgC = barcodeGen(cGenAll{ix}.idx);

% Get block representation
[outBarsCorrectLengthRescale] = barcodes_to_scaled(barcodeGen(cGenAll{ix}.idx),cGenAll{ix}.comparisonStruct, []);

consensusToPlot1 = outBarsCorrectLengthRescale(2:end,:);

% import Zeromodel.beta_cdf;
% pvalfunLeftover = @(x,l1,nuF) 1-beta_cdf(x,nuF*l1,1,0);


% pairwise compare. Need only once
allscores = zeros(size(consensusToPlot1,1),size(consensusToPlot1,1));
for i=1:size(consensusToPlot1,1)-1
    for j=i+1:size(consensusToPlot1,1)        
        overlapPx = ~isnan(consensusToPlot1(i,:)).*~isnan(consensusToPlot1(j,:));
        if sum(overlapPx) > wmin %if overlap more than a given. Could also convert to Stouffer here
            allscores(i,j) = zscore(consensusToPlot1(i,logical(overlapPx)),1)*zscore(consensusToPlot1(j,logical(overlapPx)),1)'/sum(overlapPx);
%             allscores(i,j) = norminv(1-pvalfunLeftover( allscores(i,j),sum(overlapPx),alphaNu));
        else
%             allscores(i,j) = 0.01; % tempfix - small score so consensus would still work
        end
    end
end


allscores2 = allscores'+allscores;
% square form
y = squareform(allscores2);


Z = linkage(1-y); % this links groups, but we also want to know which elements exactly are linked
% 
% numClusters = 5;
T = cluster(Z, 'cutoff', thresh,'Criterion','distance');
% [vals ids] = unique(T)
clusterBarcodes = find(T==mode(T));
% 
allScoresUp = allscores2(clusterBarcodes,clusterBarcodes);
yUp = squareform(allScoresUp);

% T = cluster(Z, "maxclust", 5);
Z = linkage(1-yUp); % this links groups, but we also want to know which elements exactly are linked


% create groups for the scores
[groups,curElts] = consensus_groups(allScoresUp,bgC(clusterBarcodes),Z);

% Now run the amplitude adjustment
[aaBarcodeIslandBlockRep,~,~,~] = consensus_cluster(consensusToPlot1(clusterBarcodes,:),Z,barcodeGen(cGenAll{ix}.idx(clusterBarcodes)),curElts,groups);

rawBarcodeIslandBlockRep = [consensusToPlot1(clusterBarcodes(groups{Z(end,1)}),:) ;consensusToPlot1(clusterBarcodes(groups{Z(end,2)}),:) ];

end

