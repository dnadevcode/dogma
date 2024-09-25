function [g] = pair_evaluation_plot(barStruct, overlapStruct,k,iy,comparisonStruct,lengthBorders, ifsave,g)
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

    if nargin < 7 
        ifsave = 0;
    end

    if iscell(barStruct)
        barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barStruct,'un',false);...
        cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end
    bBar = imresize(barStruct(k).rawBarcode(barStruct(k).rawBitmask),'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
    aBar = barStruct(iy).rawBarcode(barStruct(iy).rawBitmask);

    pA = overlapStruct(k,iy).pA ;
    pB = overlapStruct(k,iy).pB ;
    h = overlapStruct(k,iy).h;
    orr = overlapStruct(k,iy).or ;
    bS =  overlapStruct(k,iy).bestBarStretch;
       
    if nargin < 8
        f=figure;
        g=tiledlayout(2,2,'TileSpacing','compact');
    end
    nexttile(g,[1,2]);hold on

if orr==-1
    bBar = fliplr(bBar);
end

plot(-pA+1:-pA+length(bBar),zscore(bBar,1)+5,'black')
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

plot(-st+1:stop-1,zscore( bBar(pA-st+1:pA+stop-1))-6,'black')
plot(-st+1:stop-1, zscore(  aBar(pB-st+1:pB+stop-1),1)-6,'red')

% zscore(bBar(pA-st+1:pA+stop-1),1)*zscore(aBar(pB-st+1:pB+stop-1),1)'/length( bBar(pA-st+1:pA+stop-1))


    text(h+50,0,strcat('local alignment C= ', num2str( overlapStruct(k,iy).score)))
    text(stop+50,-6,strcat('full overlap C= ',num2str( overlapStruct(k,iy).fullscore)))

% [rescaleFactor, orSign]

    lgd1 = legend({strcat(['$$bar_{' num2str(k) '}$$ sF='  num2str(bS) ' or= '  num2str(orr) ]),num2str(iy)},'location','southeastoutside','Interpreter','latex')
        lgd1.Layout.Tile = 'south';

%     if orSign==-1
%         legend({num2str(ts2),strcat(['$$\bar ' num2str(bar1Name) '$$'])},'location','eastoutside','Interpreter','latex')
%     else 
%         legend({num2str(ts2),num2str(bar1Name)},'location','eastoutside')
%     end
%     stopPos = min(pos1(end),pos2(end));
%     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';

    nexttile(g,[1 2])

        markers = ['o';'s';'x';'+';'d';'v'];
    
        sets.genConsensus = 0;
    
        scaleF = 1;
    
    
    cumLengths = [0 lengthBorders]'; % plot all the theories on x axis
    
    % which theory
    posShift = cumLengths(cellfun(@(x) x.idx,comparisonStruct)');
    
    % position along the theory
    pos = cell2mat(cellfun(@(x) x.pos(1:min(end,3)),comparisonStruct,'UniformOutput',0)');
    ornt = cell2mat(cellfun(@(x) x.or(1:min(end,3)),comparisonStruct,'UniformOutput',0)');
    score = cell2mat(cellfun(@(x) x.maxcoef(1:min(end,3)),comparisonStruct,'UniformOutput',0)');
    lengthMatch = cell2mat(cellfun(@(x) x.lengthMatch,comparisonStruct,'UniformOutput',0)');

%     p3 = plot((pos+posShift)/scaleF,1:size(pos,1),'ob');
    p3(1).Marker = markers(1);
    try
        p3(2).Marker = markers(2);
    catch
    end
    try
        p3(3).Marker = markers(3);
    catch
    end
    hold on
    

%     p4.Marker = markers(1);
%         p4 = plot(pos+posShift,1:size(pos,1),'ob');

  
    if  sets.genConsensus == 1
        plot((0:100:sum(lengthBorders))/scaleF, 0.5+repmat(numBar,length(0:100:sum(lengthBorders)),1))
    end
    
    p = zeros(1,3);
    p(1) = plot(lengthBorders/scaleF,zeros(1,length(lengthBorders)),'redx');
    
    % add lines for barcode lengths
    posEnd = cell2mat(cellfun(@(x) x.pos(1)+x.lengthMatch,comparisonStruct,'UniformOutput',0)');
    for i=1:size(pos,1)
        if i==1
            colorV = 'black';
        else
           colorV = 'r';
        end
        if ornt(i)==1
            p(2) = plot(([pos(i,1) posEnd(i)]+posShift(i))/scaleF,[i i],['>',colorV,'-']);     
        else
            p(3) = plot(([pos(i,1) posEnd(i)]+posShift(i))/scaleF,[i i],['<',colorV,'-']);     
        end
    end
    
    % we also plot colorbar for best score!
    
    if scaleF == 1
        xlabel('Best position (px)','Interpreter','latex')  
    else
        xlabel('Best position (Mbp)','Interpreter','latex')  
    end
    
    legendList = {'Islands separator','Barcodes or=1','Barcodes or=2'};
    
    ylabel('Barcode nr.','Interpreter','latex')
%     if size(pos,2) == 1
        lgd2 = legend(p(find(p)), legendList(find(p)) ,'Location','southoutside','Interpreter','latex');
        lgd2.Layout.Tile = 'south';
%     else
%         legend({'$\hat C$','$C_2$','$C_3$','Theory edge seperator'},'Location','southoutside','Interpreter','latex')
%     end
    
    ylim([0,size(pos,1)+2])
   
    grid on
    grid minor
%     ax2=nexttile
%     scaledScore = score(:,1)./sqrt(lengthMatch);
%     
%     hHM =heatmap([score(:,1) scaledScore./max(scaledScore)])
%     hHM.NodeChildren(3).YDir='normal';                   % turn Y-Axis normal directionend
% %     hHM.XDisplayLabels = nan(size(hHM.XDisplayData));
% %     hHM.YDisplayLabels = nan(size(hHM.YDisplayData));
%     xlabel('score')
% %     hHM.CellLabelColor='none';
% %     linkaxes([ax1 ax2])











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

