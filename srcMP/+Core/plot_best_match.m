function [f,pos1, bar1, pos2, bar2,pcc,bar1Name,rescaleFactor, orSign] = plot_best_match(ix,stridx,mpI,LOCS,pksUniquePos,k,baridx,sF,barcodeGenGood1,h,barcodeGenGood2)
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

if nargin < 11
    barcodeGenGood2 = barcodeGenGood1; % A=B
end
% ix = 4
% In short: A - long barcode, i.e. all possible
% B - short barcode (with all possible re-scale factors)

pos = LOCS(pksUniquePos(ix)); % position on MP.
ts2 = k;

import Core.get_two_bars_alignment_params;
[bar1Name, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{ts2},mpI, pos, baridx,sF);

% print out  
% [rescaleFactor, orSign]

% maybe consider more reasonable names
b2 = barcodeGenGood1(bar1Name).rawBarcode(logical(barcodeGenGood1(bar1Name).rawBitmask));
b = barcodeGenGood2(ts2).rawBarcode(logical(barcodeGenGood2(ts2).rawBitmask));

    f=figure;nexttile([1,2]);hold on

% now we rescale b
bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
if orSign==-1
%     bRescaled = fliplr(bRescaled);
    
    bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
    pAnew = pA; % remove pA and h, this should give correct position ( check for +1 )
    plot(-pAnew+1:-pAnew+length(bRescaled),zscore(fliplr(bRescaled))+5,'black')
    pBnew = pB; % remove pA and h, this should give correct position ( check for +1 )
    plot(-pBnew+1:-pBnew+length(b2),zscore(b2)+7,'red')
%     pA
%     pB
    bar1 = fliplr(bRescaled);

%     f=figure;hold on
%     pAnew = length(bRescaled)-pA-h; % remove pA and h, this should give correct position ( check for +1 )
%     plot(-pAnew+1:-pAnew+length(bRescaled),zscore(bRescaled)+5,'black')
%     pBnew = length(b2)-pB-h; % remove pA and h, this should give correct position ( check for +1 )
%     plot(-pBnew+1:-pBnew+length(b2),zscore(fliplr(b2))+7,'red')

else
%     f=figure;hold on
    plot(-pA+1:-pA+length(bRescaled),zscore(bRescaled)+5,'black')
    plot(-pB+1:-pB+length(b2),zscore(b2)+7,'red')
    bar1 = bRescaled;


end
    
    pos1 = -pA+1:-pA+length(bRescaled);
    pos2 = -pB+1:-pB+length(b2);
    bar2 = b2;

    fP1 = find(pos1==1,1,'first');
    fP2 = find(pos2==1,1,'first');
    plot(zscore(bar1(fP1:fP1+h-1)),'black')
    plot(zscore(bar2(fP2:fP2+h-1)),'red')
    
    pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';
    pcc

    
    % now calculate PCC for total overlap
    
    [C,IA,IB] = intersect(pos1,pos2);
    pcc2 = 1/length(IA) * zscore(bar1(IA),1)*zscore(bar2(IB),1)';
    pcc2

%    bar2(IB)
%     figure;hold on
    plot(C,zscore(bar1(IA))-6,'black')
    plot(C, zscore( bar2(IB))-6,'red')
    text(h+50,0,strcat('local alignment C= ', num2str(pcc,3)))
    text(max(C)+50,-6,strcat('full overlap C= ',num2str(pcc2,3)))

% [rescaleFactor, orSign]

    legend({strcat(['$$bar_{' num2str(ts2) '}$$ '  num2str(rescaleFactor) ' '  num2str(orSign) ]),num2str(bar1Name)},'location','southeastoutside','Interpreter','latex')

%     if orSign==-1
%         legend({num2str(ts2),strcat(['$$\bar ' num2str(bar1Name) '$$'])},'location','eastoutside','Interpreter','latex')
%     else 
%         legend({num2str(ts2),num2str(bar1Name)},'location','eastoutside')
%     end
%     stopPos = min(pos1(end),pos2(end));
%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';


    saveas(f,'figs/fig3.png')


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

