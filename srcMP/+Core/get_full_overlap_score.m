function [fullscore,overlaplen,lenB, lenA, partialScore,parLen,aM,bM,aP,bP] = get_full_overlap_score(pA,pB,bestBarStretch,orr,barStruct, h)
    %get_full_overlap_score
    %
    %   Args:
    %       pA - position on A
    %       pB - position on B
    %       bestBarStretch - best bar strech
    %       orr - orientation
    %       barStruct - structure with barcode

    %   Returns:
    %       fullscore - full overlap score
    %
    %       partialScore - partial score
    %       parLen - partial length
    %       aM,bM,aP,bP - values

    if iscell(barStruct)
        barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barStruct,'un',false);...
        cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end

   % get the two barcodes (check if correct re-size method for bar1)
   aBar = imresize(barStruct(1).rawBarcode(barStruct(1).rawBitmask),'Scale' ,[1 bestBarStretch]);
   bBar = barStruct(2).rawBarcode(logical(barStruct(2).rawBitmask));

    % check if need to flip
    if orr~=1
        aBar = fliplr(aBar);
    end
    % this should be local overlap score
    %  zscore(aBar(pA:pA+h-1),1)*zscore(bBar(pB:pB+h-1),1)'/h

    % do full overlap
    lpA = length(aBar); lpB = length(bBar);

    st = min(pA,pB); % which barcode is more to the left
    stop = min(lpA-pA+1,lpB-pB+1);

    aFul = bBar(pB-st+1:pB+stop-1);
    bFul = aBar(pA-st+1:pA+stop-1);

    fullscore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
    overlaplen = length(aFul);
    
    lenB = length(bBar);
    lenA = length(aBar);

    if nargout > 6
        aM = aFul(st:st+h-1);
        bM = bFul(st:st+h-1);
    end
    
    if nargout > 4
        % now mask the common part
        aFul(st:st+h-1) = [];
        bFul(st:st+h-1) = [];
        partialScore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
        parLen = length(aFul);
    end

    if nargout > 8
        aP = aFul;
        bP = bFul;
    end

end

