function [barMat,bars,orBars,reFac] = plot_bargroup(ix,stridx,mpI,LOCS,pksUnique1,pksUniquePos,k, baridx,sF,barcodeGenGood,h,plotRes)

% here we plot a bargroup. Later this will be put into bargroup class with
% operations (such as bargroup addition/averaging)

% ix = 4
% In short: A - long barcode, i.e. all possible
% B - short barcode (with all possible re-scale factors)

barMat = cell(length(pksUniquePos),1);

b2 = barcodeGenGood(ix).rawBarcode(barcodeGenGood(ix).rawBitmask);

barMat{1}{1} = 1:length(b2);
barMat{1}{2} = b2;
barMat{1}{3} = ix;

orBars = zeros(1,length(pksUniquePos));
reFac = zeros(1,length(pksUniquePos));

if nargin>=12
    f=figure;plot(zscore(b2,1),'black'); hold on
    title('Bargroup')
end
for i=1:length(pksUniquePos)
    % check if the barcode exists
    barmatchID = find(pksUnique1{i}==ix);
    if ~isempty(barmatchID)
        pos = LOCS{i}(pksUniquePos{i}(barmatchID(1))); % position on MP.

        ts2 = k(i);


        import Core.get_two_bars_alignment_params;
        [bar1, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{ts2},mpI{i}, pos, baridx{ts2},sF);


        % maybe consider more reasonable names
        b = barcodeGenGood(ts2).rawBarcode(barcodeGenGood(ts2).rawBitmask);

        bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
        
        orBars(i) = orSign;
        reFac(i) = rescaleFactor;
        if orSign==-1
            bRescaled = fliplr(bRescaled);
    %         pA
    %         pB
        end
    
%         pAnew = length(bRescaled)-pA-h; % remove pA and h, this should give correct position ( check for +1 )
%         pBnew = length(b2)-pB-h; % remove pA and h, this should give correct position ( check for +1 )
%         plot(-pAnew+1:-pAnew+length(bRescaled),zscore(bRescaled)+5,'black')

%     plot(-pBnew+1:-pBnew+length(b2),zscore(fliplr(b2))+7,'red')

% elsee
    if nargin>=12
        plot(-pA+1+pB:-pA+pB+length(bRescaled),zscore(bRescaled)+3*i,'red')
    end
        barMat{i+1}{1} = -pA+1+pB:-pA+pB+length(bRescaled);
        barMat{i+1}{2} = bRescaled;
        barMat{i+1}{3} = ts2;

%     plot(-pB+1:-pB+length(b2),zscore(b2)+7,'red')
    else
        plot(1,1)
    end

%     end

    
end

bars = [num2str(ix) arrayfun(@(x) num2str(x),k,'un',false)];
if nargin>=12
    legend(bars,'location','southeastoutside');
    saveas(f,'figs/bargroupexample.png')
end
%% TODO: check that PCC values are correct, i.e. the same as in MP
% 
% % print out  
% [rescaleFactor, orSign]
% 
% 
% % now we rescale b
% bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
% if orSign==-1
% %     bRescaled = fliplr(bRescaled);
%     
%     bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);
% 
%     f=figure;hold on
%     pAnew = length(bRescaled)-pA-h; % remove pA and h, this should give correct position ( check for +1 )
%     plot(-pAnew+1:-pAnew+length(bRescaled),zscore(bRescaled)+5,'black')
%     pBnew = length(b2)-pB-h; % remove pA and h, this should give correct position ( check for +1 )
%     plot(-pBnew+1:-pBnew+length(b2),zscore(fliplr(b2))+7,'red')
% 
% else
%     f=figure;hold on
%     plot(-pA+1:-pA+length(bRescaled),zscore(bRescaled)+5,'black')
%     plot(-pB+1:-pB+length(b2),zscore(b2)+7,'red')
% 
% 
% end
%     
% % 
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

