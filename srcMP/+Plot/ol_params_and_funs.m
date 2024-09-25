function [aBar,bBar, subBarA, subBarB, aFul, bFul, pos] = ol_params_and_funs(barStructA, barStructB, overlapStruct)
    %
    % gens all data needed to plot in a nice format

    %   Args:
    %       barStructA - barcode A structure (source)
    %       barStructB - barcode B structure (sink)
    %       overlapStruct - results of source mapped along the sink
    %
    %   Returns:
    %   aBar - re-scaled A barcode
    %   bBar - B barocde
    %   subBarA - local subbarcode A
    %   subBarB - local subbarcode B
    %   aFul - full overlap on A
    %   bFul - full overlap on B
    %   pos - positional informaiton

    pA = overlapStruct.pA ;
    pB = overlapStruct.pB ;
    if isfield(overlapStruct,'h')
        h = overlapStruct.h;
    else
        h = overlapStruct.overlaplen; % if there is no h, only full overlap has been calculated, and we have only full overlap score
    end
    orr = overlapStruct.or ;
    bS =  overlapStruct.bestBarStretch;
     
    % use interp1 in case the imresize is done in a different way!
    if isfield(overlapStruct,'h')
        aBar = imresize(barStructA.rawBarcode(barStructA.rawBitmask),'Scale' ,[1 bS]);
        aBit = ones(1,length(aBar));

        %
        bBar = barStructB.rawBarcode(barStructB.rawBitmask);
        bBit = ones(1,length(bBar));

    else % interp1
            % barcodes
            aBar = interp1(barStructA.rawBarcode, linspace(1, length(barStructA.rawBarcode), length(barStructA.rawBarcode)*bS));
            bBar = barStructB.rawBarcode;
            % bitmasks
            aBit = barStructA.rawBitmask(round(linspace(1, length(barStructA.rawBarcode), length(barStructA.rawBarcode)*bS)));
            bBit = barStructB.rawBitmask;

            aBar(~aBit) = nan;
            bBar(~bBit) = nan;
    end

    if orr==-1 || orr==2
        aBar = fliplr(aBar);
        aBit = fliplr(aBit);
    end

    % relative positions
    pos.A = [-pA+1 -pA+length(aBar)];
    pos.B = [-pB+1 -pB+length(bBar)];

    % sub-bars, only for calculate version (includes local length parameter h)
    if isfield(overlapStruct,'h')
        subBarA = aBar(pA:pA+h-1);
        subBarB = bBar(pB:pB+h-1);
        pos.lA = [1 length(subBarA)];

        pcc = zscore(subBarA,1)*zscore(subBarB,1)'/h;
    
    else
        subBarB = [nan];
        subBarA = [nan];
    end
   
    % do full overlap
    lpA = length(aBar); lpB = length(bBar);
    
    st = min(pA,pB); % which barcode is more to the left
    stop = min(lpA-pA+1,lpB-pB+1);
    
    bFul = bBar(pB-st+1:pB+stop-1); % need to include bitmask in calculating pcc!
    aFul = aBar(pA-st+1:pA+stop-1);
    
    bfullbit = bBit(pB-st+1:pB+stop-1);
    afullbit = aBit(pA-st+1:pA+stop-1);
    pcc = zscore(aFul(logical(afullbit.*bfullbit)),1)*zscore(bFul(logical(afullbit.*bfullbit)),1)'/length(find(afullbit.*bfullbit));
  

    pos.Af = [-st+1 stop-1];

end

