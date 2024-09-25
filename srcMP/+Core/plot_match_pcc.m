function [f] = plot_match_pcc(barStruct, overlapStruct,k,iy,barStruct2,pval)
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

    if iscell(barStruct)
        barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barStruct,'un',false);...
        cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end
    
    
    if nargin < 5
        barStruct2 = barStruct;
    end
%  length(interp1(barStruct(k).rawBarcode, linspace(1, length(barStruct(k).rawBarcode), round(length(barStruct(k).rawBarcode)*overlapStruct(k,iy).bestBarStretch))))
    barSink = imresize(barStruct(k).rawBarcode,'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
    barSinkBit = imresize(barStruct(k).rawBitmask,'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
    barSink(~barSinkBit) = nan;
    try 
        barSource = imresize(barStruct(iy).rawBarcode,'Scale' ,[1 overlapStruct(k,iy).bestBarStretchSource]);
        barSourcebit = imresize(barStruct(iy).rawBitmask,'Scale' ,[1 overlapStruct(k,iy).bestBarStretchSource]);
        barSource( ~barSourcebit) = nan;
    catch
        barSource = barStruct2(iy).rawBarcode;
        barSource( ~barStruct2(iy).rawBitmask) = nan;
    end
    pA = overlapStruct(k,iy).pA ;
    pB = overlapStruct(k,iy).pB ;
    h = overlapStruct(k,iy).overlaplen;
    orr = overlapStruct(k,iy).or ;
    bS =  overlapStruct(k,iy).bestBarStretch;
       

    f=figure;nexttile([1,2]);hold on

if orr==2
    barSink = fliplr(barSink);
end
% if also source flipped
try
   if  overlapStruct(k,iy).orSource == 2 
       barSource = fliplr(barSource);
   end
catch
    
end

length(barSink)
length(barSource)
plot(-pA+1:-pA+length(barSink),nanzscore(barSink)+5,'black')
plot(-pB+1:-pB+length(barSource),nanzscore(barSource)+7,'red')

l1 = -pA+1:-pA+length(barSink); l2 = -pB+1:-pB+length(barSource);
[C,cc1,cc2] = intersect(l1,l2);

%     pos1 = -pA+1:-pA+length(bBar);
%     pos2 = -pB+1:-pB+length(aBar);
%     bar2 = b2;

%     fP1 = find(pos1==1,1,'first');
%     fP2 = find(pos2==1,1,'first');
plot(l1(cc1), nanzscore(barSink(cc1)),'black')
plot(l2(cc2), nanzscore(barSource(cc2)),'red')
nonnanpos = ~isnan(barSink(cc1)).*(~isnan(barSource(cc2)));
b1= barSink(cc1);b2 = barSource(cc2);
pcc = nanzscore(b1(find(nonnanpos)))*nanzscore(b2(find(nonnanpos)))'/sum(nonnanpos)

%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';

% 
% % do full overlap
% lpA = length(bBar); lpB = length(aBar);
% 
% st = min(pA,pB); % which barcode is more to the left
% stop = min(lpA-pA+1,lpB-pB+1);
% 
% aFul = aBar(pB-st+1:pB+stop-1);
% bFul = bBar(pA-st+1:pA+stop-1);
% 
% plot(-st+1:stop-1,nanzscore( bBar(pA-st+1:pA+stop-1))-6,'black')
% plot(-st+1:stop-1, nanzscore(  aBar(pB-st+1:pB+stop-1),1)-6,'red')

% zscore(bBar(pA-st+1:pA+stop-1),1)*zscore(aBar(pB-st+1:pB+stop-1),1)'/length( bBar(pA-st+1:pA+stop-1))


    text(max(l1(cc1))+50,0,strcat(['alignment C= ' num2str( overlapStruct(k,iy).score) ', len = ' num2str(overlapStruct(k,iy).overlaplen)]))
%     text(stop+50,-6,strcat('full overlap C= ',num2str( overlapStruct(k,iy).fullscore)))

% [rescaleFactor, orSign]
    try
        legend({strcat(['$$bar_{' num2str(k) '}$$ sF='  num2str(bS) ' or= '  num2str(orr) ]),strcat([num2str(iy) ' , sF = ' num2str(overlapStruct(k,iy).bestBarStretchSource)])},'location','southeastoutside','Interpreter','latex');     
    catch
        legend({strcat(['$$bar_{' num2str(k) '}$$ sF='  num2str(bS) ' or= '  num2str(orr) ]),num2str(iy)},'location','southeastoutside','Interpreter','latex');
    end

%     if orSign==-1
%         legend({num2str(ts2),strcat(['$$\bar ' num2str(bar1Name) '$$'])},'location','eastoutside','Interpreter','latex')
%     else 
%         legend({num2str(ts2),num2str(bar1Name)},'location','eastoutside')
%     end
%     stopPos = min(pos1(end),pos2(end));
%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';

    if nargin >=6
        title(strcat(['pval = ' num2str(pval)]))
    end
%     saveas(f,'figs/fig3.png')


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

function [scored] = nanzscore(A)

scored = (A-nanmean(A))./nanstd(A,1);
end

