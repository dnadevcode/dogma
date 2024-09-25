function [consensusToPlot1,newScore,oldScore,oldScore2] = consensus_cluster(inputData, Z, barcodeGen, curElts, groups, outBarsCorrectLengthRescale)

%   Amplitude adjustment for raw consensus using single-linkage
%   hierarchical clustering
%
%   Returns:
%       consensusToPlot1 - amplitude adjusted consensus barcode
%%
consensusOldMat = inputData;
consensusToPlot1 = inputData;

% now do the linking
newScore = zeros(1,size(Z,1));
oldScore = zeros(1,size(Z,1));
oldScore2 = zeros(1,size(Z,1));

stG = length(barcodeGen); 
for i = 1:size(Z,1)
    % barcodes
    b1 = consensusToPlot1(curElts(i,1),:);
    b2 = consensusToPlot1(curElts(i,2),:);
    % mask
    mask =logical( ~isnan(b1).*~isnan(b2));

    % means & stds
    mean1 = mean(b1(mask));
    std1 = std(b1(mask),1);
    
    mean2 = mean(b2(mask));
    std2 = std(b2(mask),1);

    % avg mean & avg std
    newMean = (mean1+mean2)/2;
    newStd = sqrt((std1^2+std2^2)/2);

    stdRatio1 = newStd/std1;
    stdRatio2 = newStd/std2;

    % now update barcodes in the clusters 1 & 3. nmean1 differ for each

    consensusToPlot1(groups{Z(i,1)},:) = (consensusToPlot1(groups{Z(i,1)},:))*stdRatio1 + newMean - mean1*stdRatio1;
    consensusToPlot1(groups{Z(i,2)},:) = (consensusToPlot1(groups{Z(i,2)},:))*stdRatio2 + newMean - mean2*stdRatio2;


    if nargin >=6 &&~isempty(outBarsCorrectLengthRescale) % comparison to theory at each step to get scores
        %    
        newVec = nanmean([consensusToPlot1(groups{Z(i,1)},:) ;consensusToPlot1(groups{Z(i,2)},:) ]);
        overlapPx = ~isnan(newVec).*~isnan(outBarsCorrectLengthRescale(1,:));
        newScore(i) = zscore(newVec(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);
            
        oldBarsVec = nanmean([consensusOldMat(groups{Z(i,1)},:) ;consensusOldMat(groups{Z(i,2)},:) ]);
    
        oldScore(i) = zscore(oldBarsVec(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);
    
        % old score z-scored
    
    
      % old score z-scored
        groups2 = [consensusOldMat(groups{Z(i,1)},:) ;consensusOldMat(groups{Z(i,2)},:) ];
        groups2 = (groups2-nanmean(groups2,2))./nanstd(groups2,1,2);
        oldBarsVec2 = nanmean(groups2);
    
        oldScore2(i) = zscore(oldBarsVec2(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);

    end

%     if nargout >=5
%         consensusRaw = nanmean([consensusOldMat(groups{Z(end,1)},:) ;consensusOldMat(groups{Z(end,2)},:) ]);
%     end
    
end




