function [f] = plot_match_simple(barStruct, overlapStruct,k,iy,ifsave)
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

 
    if isfield(overlapStruct,'pA')
        pA = overlapStruct(k,iy).pA ;
        pB = overlapStruct(k,iy).pB ;
        h = overlapStruct(k,iy).h;
        orr = overlapStruct(k,iy).or ;
        bS =  overlapStruct(k,iy).bestBarStretch;
        score = overlapStruct(k,iy).score;
        fullscore = overlapStruct(k,iy).fullscore;
    
       % resizing depends on how barcodes were defined. For PCC, we resize
        % first, and then apply bitmask. For MP?
        bBar = imresize(barStruct(k).rawBarcode(barStruct(k).rawBitmask),'Scale' ,[1 bS]);

    else
        % in this case we have slightly different structure..

        pA = overlapStruct.rezMax{iy-1}{k}.secondPos(1);
        pB = overlapStruct.rezMax{iy-1}{k}.pos(1);
        h = overlapStruct.rezMax{iy-1}{k}.lengthMatch;
        orr = overlapStruct.rezMax{iy-1}{k}.or(1);
        bS = overlapStruct.bestBarStretch{iy-1};
        score = overlapStruct.rezMax{iy-1}{k}.maxcoef(1);
        fullscore = nan; % using matlab's version, we don't calculate this within function

        lenBarTested = length(barStruct(k).rawBarcode);
        bBar = interp1(barStruct(k).rawBarcode, linspace(1,lenBarTested,lenBarTested*bS));
        bBit = barStruct(k).rawBitmask(round(linspace(1,lenBarTested,lenBarTested*bS)));
        bBar(~bBit) = nan;

    end
         aBar = barStruct(iy).rawBarcode(barStruct(iy).rawBitmask);

    f=figure;nexttile([1,2]);hold on

    if orr~=1
        bBar = fliplr(bBar);
    end
    
    plot(-pA+1:-pA+length(bBar),(bBar-mean(bBar,'omitnan' ))./std(bBar,1,'omitnan' )+5,'black')
    plot(-pB+1:-pB+length(aBar),zscore(aBar)+7,'red')

    
    %     pos1 = -pA+1:-pA+length(bBar);
    %     pos2 = -pB+1:-pB+length(aBar);
    %     bar2 = b2;
    
    %     fP1 = find(pos1==1,1,'first');
    %     fP2 = find(pos2==1,1,'first');
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
    
    plot(-st+1:stop-1,(bFul-mean(bFul,'omitnan'))./(std(bFul,1,'omitnan'  ))-6,'black')
    plot(-st+1:stop-1, zscore(  aBar(pB-st+1:pB+stop-1),1)-6,'red')
    
%     zscore(bBar(pA-st+1:pA+stop-1),1)*zscore(aBar(pB-st+1:pB+stop-1),1)'/length( bBar(pA-st+1:pA+stop-1))
    aFul = aFul(~isnan(bFul));

    bFul = bFul(~isnan(bFul));
    fullscore = zscore(bFul,1)*zscore(aFul,1)'./length(bFul)

    text(h+50,0,strcat('local alignment C= ', num2str(score)))
    text(stop+50,-6,strcat('full overlap C= ',num2str(fullscore)))

% [rescaleFactor, orSign]

    legend({strcat(['$$bar_{' num2str(k) '}$$ sF='  num2str(bS) ' or= '  num2str(orr) ]),num2str(iy)},'location','southeastoutside','Interpreter','latex')

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


% 
% if orSign==-1
%     bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
% 
%     f=figure;hold on
%     pAnew = length(bRescaled)-pA-h; % remove pA and h, this should give correct position ( check for +1 )
%     plot(-pAnew+1:-pAnew+length(bRescaled),zscore(bRescaled)+5,'black')
%     pBnew = length(b2)-pB-h; % remove pA and h, this should give correct position ( check for +1 )
%     plot(-pBnew+1:-pBnew+length(b2),zscore(fliplr(b2))+7,'red')
% 
% 
%     idx = -pAnew+1:-pAnew+length(bRescaled);
%     bar = bRescaled;
%     idx2 = -pBnew+1:-pBnew+length(b2);
%     bar2 = fliplr(b2);
% 
%     stPos= find(idx==1);
%     barCut1 = zscore(bar(stPos:stPos+h-1),1);
% 
%     stPos2= find(idx2==1);
%     barCut2 = zscore(bar2(stPos2:stPos2+h-1),1);
%     pcc = 1/length(barCut1) * barCut1*barCut2';
% 
% end


end

