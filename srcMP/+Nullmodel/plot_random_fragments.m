function [outConsensus] = plot_random_fragments(barcodeGen,goodBars,theorBar,pos,bestBarStretch,orientation)


    outConsensus = nan(length(goodBars)+1,length(theorBar));
    outConsensus(1,:) = zscore(theorBar);
    lenT = length(theorBar);
    
    for k=1:length(goodBars)
        ii = goodBars(k);
        expBar = barcodeGen{ii}.rawBarcode;
        expBit = barcodeGen{ii}.rawBitmask;

        expLen = length(expBar);

        % interpolate to the length which gave best CC value
        expBar = interp1(expBar, linspace(1,expLen,expLen*bestBarStretch(ii)));
        expBit = expBit(round(linspace(1,expLen,expLen*bestBarStretch(ii))));
        expBar(~expBit)= nan;
        
        if orientation(ii,1)==2
            expBar = fliplr(expBar);
            expBit = fliplr(expBit);
        end
        posEnd = pos(ii)+length(expBar)-1;
        numEltsOver = posEnd - lenT;
        numFirst = lenT - pos(ii)+1;
        if posEnd > lenT
            seqToplot =  (expBar-nanmean(expBar))/nanstd(expBar,1)+5;
            outConsensus(k+1,[pos(ii):lenT 1:numEltsOver]) = [seqToplot(1:numFirst) seqToplot(numFirst+1:end)] ;
        else
        outConsensus(k+1,pos(ii):pos(ii)+length(expBar)-1) = (expBar-nanmean(expBar))/nanstd(expBar,1)+5;
        end
            
    end


    isnanmat = ~isnan(outConsensus(2:end,:));
    coverage = sum(isnanmat,1);

    % f=figure,plot(find(~isnan(nanmean(outConsensus(2:end,:)))),outConsensus(2:end,~isnan(nanmean(outConsensus(2:end,:))))','.'); 
    consensus = nanmean(outConsensus(2:end,:));
    % hold on
    f=figure
    title('Synthetic plot')
    tiledlayout(3,1,'TileSpacing','compact')
    nexttile
    hold on
    plot(zscore(outConsensus(1,:))-5,'red')
    plot((consensus-nanmean(consensus))/nanstd(consensus)-5,'black')
    legend({'Synthetic theory','reference based consensus'},'location','southoutside')

    nexttile
    plot(coverage)
    legend({'Coverage'})
    nexttile
    plot(1/length(consensus)*abs(((consensus-nanmean(consensus))/nanstd(consensus)-zscore(outConsensus(1,:))).^2),'black')
    legend('RMSE')
    
    mkdir('figs');
    saveas(f,'figs/fig2.png')
%     saveas(f,'ref_comparison_2.fig')

    
    
end

