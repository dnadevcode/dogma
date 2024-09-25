function [f] = plot_match_islands(barStruct, overlapStruct,kvec,iy,ifsave)
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

    if nargin < 5 
        ifsave = 0;
    end

    if iscell(barStruct)
        barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barStruct,'un',false);...
        cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end

    aBar = barStruct(iy).rawBarcode(barStruct(iy).rawBitmask);

    f = figure;% hold on
    tiledlayout(length(kvec),2);
    % aBar always from 1

    canBeCircular = 1;

    for it=1:length(kvec)
        nexttile([1,2]);hold on
        k = kvec(it)
        bBar = imresize(barStruct(k).rawBarcode(barStruct(k).rawBitmask),'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
    
        pA = overlapStruct(k,iy).pA ;
        pB = overlapStruct(k,iy).pB ;
        h = overlapStruct(k,iy).h;
        orr = overlapStruct(k,iy).or ;
        bS =  overlapStruct(k,iy).bestBarStretch;
       
        %

        if orr==-1
            bBar = fliplr(bBar);
        end

        plot(-pB+1:-pB+length(aBar),zscore(aBar)+7,'red')

        if  canBeCircular ==1 && -pA+length(bBar)>-pB+length(aBar)
            rightPx = -pA+1:-pB+length(aBar);
            leftPx = -pB+length(aBar)+1:-pA+length(bBar);
            plot(rightPx,zscore(bBar(1:length(rightPx)))+5,'black')
            plot(-pB+1:-pB+length(leftPx),zscore(bBar(length(rightPx)+1:end))+5,'black')

        else
            plot([-pA+1:-pA+length(bBar)],zscore(bBar,1)+5,'black')
        end

        plot(zscore(bBar(pA:pA+h-1),1),'black')
        plot(zscore(aBar(pB:pB+h-1),1),'red')
        
        pcc = zscore(aBar(pB:pB+h-1),1)*zscore(bBar(pA:pA+h-1),1)'/h

%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';

    
    % do full overlap
    lpA = length(bBar); lpB = length(aBar);
    
    st = min(pA,pB); % which barcode is more to the left
    stop = min(lpA-pA+1,lpB-pB+1);
    
    aFul = aBar(pB-st+1:pB+stop-1);
    bFul = bBar(pA-st+1:pA+stop-1);
    
    plot(-st+1:stop-1,zscore( bBar(pA-st+1:pA+stop-1))-6,'black')
    plot(-st+1:stop-1, zscore(  aBar(pB-st+1:pB+stop-1),1)-6,'red')
    if  canBeCircular ==1 && -pA+length(bBar)>-pB+length(aBar)
        plot(-pB+1:-pB+length(leftPx),zscore(bBar(length(rightPx)+1:end))-6,'black')
        plot(-pB+1:-pB+length(leftPx),zscore(aBar(1:length(-pB+1:-pB+length(leftPx))))-6,'red')

    end

% zscore(bBar(pA-st+1:pA+stop-1),1)*zscore(aBar(pB-st+1:pB+stop-1),1)'/length( bBar(pA-st+1:pA+stop-1))

    
    text(h+50,0,strcat('local alignment C= ', num2str( overlapStruct(k,iy).score)))
    text(stop+50,-6,strcat('full overlap C= ',num2str( overlapStruct(k,iy).fullscore)))

% [rescaleFactor, orSign]

%    lgd= legend({strcat(['$$bar_{' num2str(k) '}$$ sF='  num2str(bS) ' or= '  num2str(orr) ]),num2str(iy)},'location','southeastoutside','Interpreter','latex')
%     lgd.Location = 'southoutside';
    end
%     if orSign==-1
%         legend({num2str(ts2),strcat(['$$\bar ' num2str(bar1Name) '$$'])},'location','eastoutside','Interpreter','latex')
%     else 
%         legend({num2str(ts2),num2str(bar1Name)},'location','eastoutside')
%     end
%     stopPos = min(pos1(end),pos2(end));
%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';

if ifsave
    [~,~] = mkdir('figs');
    saveas(f,'figs/fig3.png')
end


end

