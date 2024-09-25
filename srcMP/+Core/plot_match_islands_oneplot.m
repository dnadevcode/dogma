function [f,globExp] = plot_match_islands_oneplot(barStruct,cGenAll, overlapStruct,kvec,iy,ifsave,f)
    %   plot_best_match
    %   Args:
    %       ix,stridx,mpI,LOCS,pksUniquePos,k,baridx,sF,barcodeGenGood1,h,barcodeGenGood2
    %
    %   ix - which max to select
    %   stridx - re-scaling and orientation factors
    %   mpI - matrix profile index
    %   LOCS -  
    %   pksUniquePos - positions in the LOCS vector of unique barcodes
    %   k - 
    %   baridx - barcode indexe in the concatenated time series
    %   sF - re-scaling factors
    %   barcodeGenGood1 - barcodes which are concatenated into long
    %   time-series when running local comparison
    %   barcodeGenGood2 - barcodes that are re-scaled and oriented when
    %   running local comparison
    %   h - overlap window size

    %   Returns:
    %       f,pos1, bar1, pos2, bar2,pcc,bar1Name
    %

    %   Example

    % Possible case: exp longer than theory

    if nargin < 6 
        ifsave = 0;
    end

    if iscell(barStruct)
        barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barStruct,'un',false);...
        cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end

    aBar = barStruct(iy).rawBarcode(barStruct(iy).rawBitmask);
    if nargin < 7
        f = figure;hold on
            tiledlayout(2,1);

    end
    nexttile;hold on
%     tiledlayout(length(kvec),2);
    % aBar always from 1
    plot(zscore(aBar)+7,'red')
    % -pB+1:-pB+length(aBar)
    canBeCircular = 1;

%     globExp = [] % if for theory
    globExp = zeros(sum(cellfun(@(x) length(x.idx),cGenAll)),length(aBar));
    itBar = 1;



    for it=1:length(kvec)
%         nexttile([1,2]);hold on
        k = kvec(it)
        bBar = imresize(barStruct(k).rawBarcode(barStruct(k).rawBitmask),'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
    
        pA = overlapStruct(k,iy).pA ;
        pB = overlapStruct(k,iy).pB ;
        h = overlapStruct(k,iy).h;
        orr = overlapStruct(k,iy).or ;
        bS =  overlapStruct(k,iy).bestBarStretch;
       
        %
        shiftCur =  -pB+1;

        if orr==-1
            bBar = fliplr(bBar);
        end


        if  canBeCircular == 1 && (-pA+length(bBar)>-pB+length(aBar))
            rightPx = -pA+1:-pB+length(aBar);
            leftPx = -pB+length(aBar)+1:-pA+length(bBar);
            plot(rightPx-shiftCur,zscore(bBar(1:length(rightPx)))+5+it*2,'black')
            plot([-pB+1:-pB+length(leftPx)]-shiftCur+1,zscore(bBar(length(rightPx)+1:end))+5+it*2,'black');
            globPx = [[rightPx-shiftCur] [-pB+1:-pB+length(leftPx)]-shiftCur+1];
        else
            if canBeCircular == 1 &&  -pA < -pB
                % doesn't fit on the lhs, move to rhs
                leftPx = [-pB+1:-pA+length(bBar)];

                rightPx = [-pB+length(aBar)-(-pB+1 - (-pA+1))+1:-pB+length(aBar)];
                
%                 [length(aBar)-(pA-pB+1):length(aBar)];
                
                plot(rightPx-shiftCur,zscore(bBar(1:length(rightPx)))+5+it*2,'black')
                plot([-pB+1:-pB+length(leftPx)]-shiftCur+1,zscore(bBar(length(rightPx)+1:end))+5+it*2,'black');

                globPx = [[rightPx-shiftCur] [-pB+1:-pB+length(leftPx)]-shiftCur+1];

%                 canBeCircular
            else
                plot([-pA+1:-pA+length(bBar)]-shiftCur,zscore(bBar,1)+5+it*2,'black')
                globPx = [-pA+1:-pA+length(bBar)]-shiftCur;
            end
        end

        % create coverage map from overlaps
        allPosLeft = cellfun(@(x) x.pos, cGenAll{it}.comparisonStruct);

        for kk=1:length(cGenAll{it}.comparisonStruct)
            % here place the barcode wrt global
            curPos = round((cGenAll{it}.comparisonStruct{kk}.pos-min(allPosLeft)+1)*overlapStruct(k,iy).bestBarStretch);
            globExp(itBar,globPx(max(1,curPos):min(end,curPos+round(cGenAll{it}.comparisonStruct{kk}.lengthMatch*overlapStruct(k,iy).bestBarStretch)-1))) = 1;
            itBar = itBar + 1;

        end

        plot([1:h]-shiftCur,zscore(bBar(pA:pA+h-1),1),'black')
        plot([1:h]-shiftCur, zscore(aBar(pB:pB+h-1),1),'red')
        
        pcc = zscore(aBar(pB:pB+h-1),1)*zscore(bBar(pA:pA+h-1),1)'/h

%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';

    
    % do full overlap
    lpA = length(bBar); lpB = length(aBar);
    
    st = min(pA,pB); % which barcode is more to the left
    stop = min(lpA-pA+1,lpB-pB+1);
    
    aFul = aBar(pB-st+1:pB+stop-1);
    bFul = bBar(pA-st+1:pA+stop-1);
    
    plot([-st+1:stop-1]-shiftCur,zscore( bBar(pA-st+1:pA+stop-1))-6-it*2,'black')
    plot([-st+1:stop-1]-shiftCur, zscore(  aBar(pB-st+1:pB+stop-1),1)-6-it*2,'red')
    if  canBeCircular ==1 && -pA+length(bBar)>-pB+length(aBar)
        plot([-pB+1:-pB+length(leftPx)]-shiftCur,zscore(bBar(length(rightPx)+1:end))-6-it*2,'black')
        plot([-pB+1:-pB+length(leftPx)]-shiftCur,zscore(aBar(1:length(-pB+1:-pB+length(leftPx))))-6-it*2,'red')
    end

    end
    title('A) Visual comparison')
    nexttile
    plot(sum(globExp))
    hold on
    title('B) Coverage')
    % coverage

if ifsave
    [~,~] = mkdir('figs');
    saveas(f,'figs/fig3.png')
end


end

