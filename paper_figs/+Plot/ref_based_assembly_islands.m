function [f] = ref_based_assembly_islands(barcodeGen,barIslands,synthStr,theoryStruct)

f=figure('Position', [10 10 900 800])
tiledlayout(2*length(barIslands),1,'TileSpacing','compact');


for idx = 1:length(barIslands)
    bars = barcodeGen(barIslands{idx}); 
    barx = barIslands{idx};
    goodBars = 1:length(barx);

    comparisonStruct = synthStr(barx);
% [outConsensus2, coverage2, pval3,f] = gen_reference_based_assembly(bars,synthStr{ixtest}(barx),theoryStruct{ixtest},'test11',inf);

    % convert the reference coefs to p-vals, etc.
    orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');
 
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
        expBar = bars{ii}.rawBarcode;
        expBit = bars{ii}.rawBitmask;

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

nexttile
hold on
title(['Consensus vs. theory for barcode island ', num2str(idx)],'Interpreter','latex');
plot(zscore(outConsensus(1,:))-5,'red')
plot((consensus-nanmean(consensus))/nanstd(consensus)-5,'black')
lgd1= legend({'In-silico predicted barcode','reference based consensus'});
lgd1.Location = 'eastoutside';
nexttile
plot(coverage)
lgd2 = legend({'Coverage'});
lgd2.Location = 'eastoutside';

% nexttile
% plot(1/length(consensus)*abs(((consensus-nanmean(consensus))/nanstd(consensus)-zscore(outConsensus(1,:))).^2),'black')
% lgd3 =legend('RMSE');
% lgd3.Location = 'eastoutside';



end

