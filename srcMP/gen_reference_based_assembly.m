function [outConsensus,coverage, pval,consensus,f] = gen_reference_based_assembly(barcodeGen,comparisonStruct,theoryStruct,timestamp,pvalThresh)

    if nargin < 5
        pvalThresh = 0.01; % thresh for assigning barcodes
    end 

%     MIN_OVERLAP_PIXELS = 500; %
    % convert the reference coefs to p-vals, etc.
    orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');
    maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',0)');

%     len2 = 
    
    idxthry = cell2mat(cellfun(@(x) x.idx,comparisonStruct,'UniformOutput',0)');
    lengthsThry = arrayfun(@(x) theoryStruct{x}.length,idxthry);
    
    lengthsbar = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    
    % pvalue params. should change based on psf
    ccthresh = 0.3;
    nuF = 0.13;
    import Zeromodel.beta_ev_cdf;
    pval = nan(1,length(barcodeGen));
    for i=1:length(barcodeGen)

        if comparisonStruct{i}.maxcoef(1) > ccthresh
            pval(i) = 1-beta_ev_cdf(comparisonStruct{i}.maxcoef(1), nuF*lengthsbar(i), 1, 2*max(lengthsbar(i),lengthsThry(i)),1);
        end
    end

% %     for i=1:length(barcodeGen)
% %         min_overlap = lengthsbar(i); %min(lengthsThry(i),lengthsbar(i));
% %         a = 0.13*min_overlap;
% %         n = 2*(lengthsThry(i)-min_overlap);
% %         pval(i) = 1-Zeromodel.beta_ev_cdf(maxcoef(i), a, 1, n, 1);
% % 
% % %         0.004*(len2-MIN_OVERLAP_PIXELS+1)*(len1-MIN_OVERLAP_PIXELS+1);
% %     end
    
    goodBars = find(pval< pvalThresh);
    
    % load theory file. just single one
    try
        theorBar = theoryStruct{comparisonStruct{1}.idx}.rawBarcode;
    catch
        fileID = fopen(theoryStruct{comparisonStruct{1}.idx}.filename,'r');
        formatSpec = '%f';
        theorBar = fscanf(fileID,formatSpec);
        fclose(fileID);
    end
    %% plot multi theories with some nan's in between
    outConsensus = nan(length(goodBars)+1,length(theorBar));
    outConsensus(1,:) = zscore(theorBar);
    lenT = length(theorBar);
    
    for k=1:length(goodBars)
        ii = goodBars(k);
        expBar = barcodeGen{ii}.rawBarcode;
        expBit = barcodeGen{ii}.rawBitmask;

        expLen = length(expBar);

        % interpolate to the length which gave best CC value
        expBar = interp1(expBar, linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch));
        expBit = expBit(round(linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch)));
        expBar(~expBit)= nan;
        
        if orientation(ii,1)~=1
            expBar = fliplr(expBar);
            expBit = fliplr(expBit);
        end
    
        if pos(ii) < 1
            pos(ii) = pos(ii)+lenT;
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
%%
isnanmat = ~isnan(outConsensus(2:end,:));
coverage = sum(isnanmat,1);

% f=figure,plot(find(~isnan(nanmean(outConsensus(2:end,:)))),outConsensus(2:end,~isnan(nanmean(outConsensus(2:end,:))))','.'); 
consensus = nanmean(outConsensus(2:end,:));
consensus(:,coverage<=5) = nan;

% PCC score for this
pcc = 1/sum(~isnan(consensus)) * zscore(consensus(~isnan(consensus)))*zscore(outConsensus(1,~isnan(consensus)))';

% hold on
f=figure;
tiledlayout(3,1,'TileSpacing','compact');
nexttile
hold on
plot(zscore(outConsensus(1,:))-5,'red')
plot((consensus-nanmean(consensus))/nanstd(consensus)-5,'black')
lgd1= legend({'In-silico predicted barcode','reference based consensus'});
lgd1.Location = 'eastoutside';
nexttile
plot(coverage)
lgd2 = legend({'Coverage'});
lgd2.Location = 'eastoutside';

nexttile
plot(1/length(consensus)*abs(((consensus-nanmean(consensus))/nanstd(consensus)-zscore(outConsensus(1,:))).^2),'black')
lgd3 =legend('RMSE');
lgd3.Location = 'eastoutside';

try
    saveas(f,fullfile(timestamp,'ref_comparison_2.png'))
    saveas(f,fullfile(timestamp,'ref_comparison_2.fig'))
catch
end
%     saveas(f,'figs/fig4.png')



    
    
    % now plot.
    
    

end

