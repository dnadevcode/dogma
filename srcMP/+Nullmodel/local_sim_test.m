function [curStd,pairs,scrs] = local_sim_test(pscores,barcodeGen, overlapStruct,MIN_OVERLAP_PIXELS,N,NN,timestamp,overlapStructMP)
%%
% MIN_OVERLAP_PIXELS = 150;
% 
% N = 150;
% NN = 2000;

if nargin < 8
    calcMP = 1;
else
    calcMP = 0;
end

% PCC_OVERLAP = reshape([overlapStruct1.sc], size(overlapStruct1,1),size(overlapStruct1,2));
PCC_OVERLAP = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2));

idxPair = 1;
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
% normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(filtM(:),'desc','MissingPlacement','last');

% [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

NN = sum(~isnan(sortedVals));
pairs = zeros(2,NN);
scrs = cell(1,NN);
curStd = zeros(1,NN);
for idxPair = 1:NN;
    idxPair
    
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
%     if calcMP
% %         [overlapStruct2] = calc_overlap_mp(barcodeGen([xId,yId]),overlapStruct(xId,yId).bestBarStretch, MIN_OVERLAP_PIXELS,timestamp);
%     else
%         overlapStruct2 = overlapStructMP(xId,yId);
%     end
    % get the allposB and allposA for overlapping positions
    curPts = overlapStructMP(xId,yId).fulloverlapPosARescaled-overlapStructMP(xId,yId).fulloverlapPosRoot;
%     curPts =  overlapStruct2(1,2).allposB(1:min(end,N))-overlapStruct2(1,2).allposA(1:min(end,N));
%     curStd(idxPair)  = std( overlapStruct2(1,2).allposB(1:min(end,N))-overlapStruct2(1,2).allposA(1:min(end,N)));
    curStd(idxPair) = std(curPts);
    pairs(:,idxPair) = [xId;yId];
   scrs{idxPair} = curPts;
end

% thrStd = 50;
% f=figure
% plot(cumsum(curStd>thrStd))

% 
% barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
%     cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
% 
% import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, overlapStructMP,xId,yId);
% 
% figure,
% tiledlayout(3,1)
% nexttile
% plot(overlapStructMP(xId,yId).fulloverlapor)
% nexttile
% plot(abs(overlapStructMP(xId,yId).localSFoverlap))
% nexttile
% plot( overlapStructMP(xId,yId).fulloverlapPosA-overlapStructMP(xId,yId).fulloverlapPosRoot)
% hold on 
% plot(overlapStructMP(xId,yId).fulloverlapPosARescaled-overlapStructMP(xId,yId).fulloverlapPosRoot)


end


