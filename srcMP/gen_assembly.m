function [outConsensus,coverage,consensus,islandElts,islandPx,f] = gen_assembly(barcodeGen, comparisonStruct, lengthsThry, timestamp, minCoverage,plot)

    if nargin < 5
        minCoverage = 5; % minimum coverage
    end
%     MIN_OVERLAP_PIXELS = 500; %
    % convert the reference coefs to p-vals, etc.
    orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');
    maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',0)');
    strF = cell2mat(cellfun(@(x) x.bestBarStretch,comparisonStruct,'UniformOutput',0)');

    idxthry = cell2mat(cellfun(@(x) x.idx,comparisonStruct,'UniformOutput',0)');
    lengthsbar = cellfun(@(x) length(x.rawBarcode),barcodeGen);

    outConsensus = cell(1,length(unique(idxthry)));
    coverage = cell(1,length(unique(idxthry)));
    consensus = cell(1,length(unique(idxthry)));
    islandElts = cell(1,length(unique(idxthry)));
    islandPx = cell(1,length(unique(idxthry)));
    %%
    for ix = 1:length(unique(idxthry))
        curIx = find(idxthry==ix);
        posCur = pos(curIx);
%         maxcoefCur = maxcoef(curIx);
        strFCur = strF(curIx);
        orientationCur = orientation(curIx);
        lengthsbarCur = lengthsbar(curIx);

        posCur = posCur-min(posCur)+1;
   %% plot multi theories with some nan's in between
        maxLen = max(posCur+round(lengthsbarCur.*strFCur));
        outConsensus{ix} = nan(length(lengthsbarCur),maxLen);
%     outConsensus(1,:) = zscore(theorBar);
%     lenT = length(theorBar);
    
        for k=1:length(lengthsbarCur)
            expBar = barcodeGen{curIx(k)}.rawBarcode;
            expBit = barcodeGen{curIx(k)}.rawBitmask;

            % interpolate to the length which gave best CC value
            expBar = imresize(expBar,'Scale' ,[1 strFCur(k)]);
            expBit = imresize(expBit,'Scale' ,[1 strFCur(k)]);
            expBar(~expBit)= nan;

            if orientationCur(k,1)==2
                expBar = fliplr(expBar);
                expBit = fliplr(expBit);
            end
            outConsensus{ix}(k,posCur(k):posCur(k)+length(expBar)-1) = (expBar-nanmean(expBar))/nanstd(expBar,1)+5;
        end
% outConsensus,coverage,consensus,islandElts,islandPx
    %%
    isnanmat = ~isnan(outConsensus{ix}(1:end,:));
    coverage{ix} = sum(isnanmat,1);

    consensus{ix} = nanmean(outConsensus{ix}(1:end,:));
    consensus{ix}(:,coverage{ix} <= minCoverage) = nan;

    islands = bwconncomp(~isnan(consensus{ix}));

    idxs = cellfun(@(x) find(~isnan(nanmean(outConsensus{ix}(:,x)'))),islands.PixelIdxList,'un',false);
    islandElts{ix} = curIx(idxs{1}) ;
    islandPx{ix} = cellfun(@(x) x ,islands.PixelIdxList,'un',false);

    % bars convering a specific island

    % PCC score for this
%     pcc = 1/sum(~isnan(consensus)) * zscore(consensus(~isnan(consensus)))*zscore(outConsensus{ix}(1,~isnan(consensus)))';

    end

    if nargin >= 6
       % hold on
        f=figure;
        tiledlayout(2,2,'TileSpacing','compact')
        ax1=nexttile
        hold on
        % plot(zscore(outConsensus(1,:))-5,'red')
        for ix=1:length(consensus)
            for ii=1:length(islandPx{ix})
                plot(islandPx{ix}{ii},(consensus{ix}(islandPx{ix}{ii})-nanmean(consensus{ix}(islandPx{ix}{ii})))/nanstd(consensus{ix}(islandPx{ix}{ii}))-5*ix)
            end
        end
        title('Consensus')
    %     legend({'Barcode islands'},'location','southoutside')

        ax2=nexttile
        title('Coverage')
        hold on
        for ix=1:length(consensus)
            plot(coverage{ix})
        end
    %     legend({'Coverage'},'location','best')
        ax3=nexttile
        hold on
        title('Standard deviation')
        for ix=1:length(consensus)
            plot(nanstd(outConsensus{ix})+ix)
        end
    %     legend('STD')
        linkaxes([ax1 ax2 ax3],'x')
    %     nexttile
        legVals = arrayfun(@(x) strcat('Consensus\_barcode\_',num2str(x)),1:length(outConsensus),'un',false);
        lgd = legend( legVals);
        lgd.Layout.Tile = 4;
        try
            saveas(f,fullfile(timestamp,'ref_comparison_2.png'))
            saveas(f,fullfile(timestamp,'ref_comparison_2.fig'))
        catch
        end
    end

    % f=figure,plot(find(~isnan(nanmean(outConsensus(2:end,:)))),outConsensus(2:end,~isnan(nanmean(outConsensus(2:end,:))))','.'); 


 
    %     saveas(f,'figs/fig4.png')




%     len2 = 
    
%     lengthsThry = arrayfun(@(x) theoryStruct{x}.length,idxthry);
    
    
%     ccthresh = 0.3;
%     nuF = 0.03;
%     import Zeromodel.beta_ev_cdf; % do we need here? already best barcodes by sorted values..
%     pval = nan(1,length(barcodeGen));
%     for i=1:length(barcodeGen)
% 
%         if comparisonStruct{i}.maxcoef > ccthresh % this is for overlap, so change..
%             pval(i) = 1-beta_ev_cdf(comparisonStruct{i}.maxcoef, nuF*lengthsbar(i), 1, 2*max(lengthsbar(i),lengthsThry(idxthry(i))),1);
%         end
%     end

% %     for i=1:length(barcodeGen)
% %         min_overlap = lengthsbar(i); %min(lengthsThry(i),lengthsbar(i));
% %         a = 0.13*min_overlap;
% %         n = 2*(lengthsThry(i)-min_overlap);
% %         pval(i) = 1-Zeromodel.beta_ev_cdf(maxcoef(i), a, 1, n, 1);
% % 
% % %         0.004*(len2-MIN_OVERLAP_PIXELS+1)*(len1-MIN_OVERLAP_PIXELS+1);
% %     end
    
%     goodBars = find(pval< 0.01);
    
    % load theory file. just single one
%     fileID = fopen(theoryStruct{comparisonStruct{1}.idx}.filename,'r');
%     formatSpec = '%f';
%     theorBar = fscanf(fileID,formatSpec);
%     fclose(fileID);
    


    
    
    % now plot.
    
    

end

