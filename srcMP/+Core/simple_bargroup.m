function [orBars,reFac,pA,pB,barsFinal] = simple_bargroup(barStruct, overlapStruct,root,pthresh,toplot)
%   Creates simple bargroup for the specified input.

% todo: for this keep only the edges between the considered barcodes. Can
% create overlap graph.

%   Args:
%     PCC_OVERLAP,PCC_MP,len1,len2,lenOverlap

%   Returns:
%       

PCC_OVERLAP = reshape([overlapStruct.fullscore], size(overlapStruct,1),size(overlapStruct,2));
lenRoot = reshape([overlapStruct.lenB], size(overlapStruct,1),size(overlapStruct,2));
lenA = reshape([overlapStruct.lenA], size(overlapStruct,1),size(overlapStruct,2));
overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));

% index of ROOT barcode 
pccVals = PCC_OVERLAP(:,root); % columns are non re-scaled

% overlapStruct
%%
% bar1 = 20;
% figure,plot(PCC_MP(bar1,:))
% figure,plot(PCC_OVERLAP(bar1,:))

idxToCheck = find(pccVals>0.5); % check only id's passing initial thresh

pval = nan(1,length(pccVals));
import Zeromodel.beta_ev_pdf;
for j=idxToCheck'
        a = 0.08*overlaplen(j,root); % based on pval_test tests
        n = 2*(lenRoot(j, root)+lenA(j,root)-2*overlaplen(j,root));
        pval(j) = 1-Zeromodel.beta_ev_cdf( PCC_OVERLAP(j,root), a, 1, n, 1);
end

numGoodMatch = pval < pthresh;

goodvals = pval(pval < pthresh);
peaksToTry = find(numGoodMatch);

orBars = [overlapStruct(peaksToTry,root).or];
reFac = [overlapStruct(peaksToTry,root).bestBarStretch];
pA = [overlapStruct(peaksToTry,root).pA];
pB = [overlapStruct(peaksToTry,root).pB];
barsFinal=[root peaksToTry];

if ~isempty(toplot)
    fig = figure; 
    tiledlayout(1,2);
    nexttile([1 2])
    hold on
    plot(zscore(barStruct(root).rawBarcode(barStruct(root).rawBitmask)),'black')
    for j=1:length(peaksToTry)
        rescaled = imresize(barStruct(peaksToTry(j)).rawBarcode(barStruct(peaksToTry(j)).rawBitmask),'Scale' ,[1 reFac(j)]);
        if orBars(j)==-1
            rescaled = fliplr(rescaled);
        end
        plot(pB(j)-pA(j)+1:pB(j)-pA(j)+length(rescaled),zscore(rescaled)+3*j,'black')
    end

    bars = cell(1,length(peaksToTry)+1);
    bars{1} = strcat('root ',num2str(root));
%     [num2str(root) arrayfun(@(x) num2str(x),[peaksToTry],'un',false)];
    for j=1:length(peaksToTry)
         bars{j+1} = strcat([num2str(peaksToTry(j)) '; or=' num2str( orBars(j))...
             '; rf=' num2str(reFac(j)) '; pcc=' num2str(pccVals(peaksToTry(j)),'%3.2f') '; pval=' num2str(goodvals(j),'%4.2e')])
    end
        legend(fliplr(bars),'location','southoutside','Interpreter','latex');
        saveas(fig,'figs/fig6.png')

%      for j=1:length(peaksToTry)+1
%          plot( )

%     if nargin>=12
%         saveas(f,'figs/bargroupexample.png')
%     end

end
% 
%     % bargroup..
%     import Core.plot_bargroup;
%     if pl==1
%         [barMat,bars,orBars,reFac] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,pl);
%     else
%         [barMat,bars,orBars,reFac] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
%     end        
%     % instead create matrix like in gen reference based assembly
% %     figure;hold on;
% %     for i=1:length(barMat)
% %         plot(  barMat{i}{1},zscore(  barMat{i}{2} )+3*i,'red');
% %     end

end

