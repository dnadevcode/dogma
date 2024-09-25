
% given theory barcode, we investigate variation (over w=300)
w = 300;

thrMean = arrayfun(@(x) mean(theoryStruct{1}{1}.rawBarcode(x:x+w-1)),1:length(theoryStruct{1}{1}.rawBarcode)-w+1);
thrStd = arrayfun(@(x) std(theoryStruct{1}{1}.rawBarcode(x:x+w-1)),1:length(theoryStruct{1}{1}.rawBarcode)-w+1);

figure;tiledlayout(2,1);
nexttile
plot(thrMean)
title('Moving mean theory')
nexttile
plot(thrStd)
title('Moving std theory')
xlabel('px')


% w = 300;
% 
% thrMean = arrayfun(@(x) mean(theoryStruct.rawBarcode(x:x+w-1)),1:length(theoryStruct.rawBarcode)-w+1);
% thrStd = arrayfun(@(x) std(theoryStruct.rawBarcode(x:x+w-1)),1:length(theoryStruct.rawBarcode)-w+1);
% 

% Create properly z-scored consensus (no bias from different AT/GC regions)
% import Core.create_consensus;
% [consensusBar] = create_consensus(barcodeIslandsData{1}, barcodeGen(cGenAll{1}.idx));


% Reference based assembly - this is just makes everything zero mean std 1
[outConsensus2, coverage2, pval3,f] = gen_reference_based_assembly(barcodeGen,synthStr{1},theoryStruct{1},'test11',inf);


% Tobias/Albertas idea
[outBarsCorrectLengthRescale] = barcodes_to_scaled(barcodeGen,synthStr{1},theoryStruct{1});

expBars = outBarsCorrectLengthRescale(2:end,:);

% sort the data based on starting position
barsNonCrossingZero = isnan(expBars(:,1));

 posidx = zeros(1,size(expBars,1));
posidx(barsNonCrossingZero) = arrayfun(@(x) find(~isnan(expBars(x,:)),1,'first'),find(barsNonCrossingZero));
posidx(~barsNonCrossingZero)= arrayfun(@(x) find(isnan(expBars(x,:)),1,'last'),find(~barsNonCrossingZero));
[pos, idxv] = sort(posidx,'desc');
% pos(barsNonCrossingZero) = idxv;
% pos(barsNonCrossingZero) = idxv2;


% [pos2,idxv2] = sort(arrayfun(@(x) find(isnan(expBars(x,:)),1,'first')+size(expBars,2) +(find(isnan(expBars(x,:)),1,'last')-find(~isnan(expBars(x,:)),1,'first'))/2,find(~barsNonCrossingZero)));


% [pos,idxv] = sort(arrayfun(@(x) find(~isnan(expBars(x,:)),1,'first') +(find(~isnan(expBars(x,:)),1,'last')-find(~isnan(expBars(x,:)),1,'first'))/2,1:size(expBars,1)));
% 
consensusToPlot1 = expBars(idxv,:);
figure,imagesc(consensusToPlot1);colormap(gray)
% figure,imagesc(expBars)

% pairwise compare
allscores = zeros(size(consensusToPlot1,1),size(consensusToPlot1,1));
for i=1:size(consensusToPlot1,1)-1
    for j=i+1:size(consensusToPlot1,1)        
        overlapPx = ~isnan(consensusToPlot1(i,:)).*~isnan(consensusToPlot1(j,:));
        if sum(overlapPx) > 300 %if overlap more than a given. Could also convert to Stouffer here
            allscores(i,j) = zscore(consensusToPlot1(i,logical(overlapPx)),1)*zscore(consensusToPlot1(j,logical(overlapPx)),1)'/sum(overlapPx);
        end
    end
end

allscores2 = allscores'+allscores;
% square form
y = squareform(allscores2);


Z = linkage(1-y); % this links groups, but we also want to know which elements exactly are linked

allscores(allscores==0) = nan;
[sortedVals, sortedIds] = sort(allscores(:),'desc','MissingPlacement','last');
[curSinkVec,curSourceVec] =  arrayfun(@(x) ind2sub(size(allscores),sortedIds(x)),1:numel(sortedVals));
curVec = [curSinkVec;curSourceVec];

% get the first merging element of each group
groups = [num2cell(1:length(barcodeGen)) cell(1,size(Z,1))];
curElts = zeros(size(Z,1),2);
stG = length(barcodeGen);
st = 1;
for i = 1:size(Z,1)
    groups{stG+i} = [groups{Z(i,1)} groups{Z(i,2)}];
    % we find first element in curVec that has two elements from this group
    for j=st:size(curVec,2)
        s1 = sum(ismember(curVec(:,j),groups{Z(i,1)}))==1;
        s2 = sum(ismember(curVec(:,j),groups{Z(i,2)}))==1;
        if s1&&s2
            if ismember(curVec(1,j),groups{Z(i,1)})   
                curElts(i,:) = curVec(:,j);
            else
                curElts(i,:) = curVec([2 1],j);
            end
            st = st+1;
            break;
        end
    end
end


figure % can plot the dendogram
dendrogram(Z); title('Barcode dendogram')

%%
consensusOldMat =  expBars(idxv,:);
consensusToPlot1 = expBars(idxv,:);

% now do the linking
newScore = zeros(1,size(Z,1));
oldScore = zeros(1,size(Z,1));
oldScore2 = zeros(1,size(Z,1));

stG = length(barcodeGen); 
for i = 1:size(Z,1)
    % barcodes
    b1 = consensusToPlot1(curElts(i,1),:);
    b2 = consensusToPlot1(curElts(i,2),:);
    % mask
    mask =logical( ~isnan(b1).*~isnan(b2));

    % means & stds
    mean1 = mean(b1(mask));
    std1 = std(b1(mask),1);
    
    mean2 = mean(b2(mask));
    std2 = std(b2(mask),1);

    meano1 = mean(b1(~mask),'omitnan');
    stdo1= std(b1(~mask),1,'omitnan');
    
    meano2 = mean(b2(~mask),'omitnan');
    stdo2 = std(b2(~mask),1,'omitnan');


    % means total
    nmean1 = mean(consensusToPlot1(groups{Z(i,1)},mask),2,'omitnan');
    nmean2 = mean(consensusToPlot1(groups{Z(i,2)},mask),2,'omitnan');
    nmeano1 = mean(consensusToPlot1(groups{Z(i,1)},~mask),2,'omitnan');
    nmeano2 = mean(consensusToPlot1(groups{Z(i,2)},~mask),2,'omitnan');
    nstd1 = std(b1,1,'omitnan');
    nstd2 = std(b2,1,'omitnan');

    % avg mean & avg std
    newMean = (mean1+mean2)/2;
    newStd = sqrt((std1^2+std2^2)/2);

    meanRatio1 = newMean/mean1;%*nmean1/mean1;
    stdRatio1 = newStd/std1;


    meanRatio2 = newMean/mean2;%*nmean2/mean2;
    stdRatio2 = newStd/std2;

    % now update barcodes in the clusters 1 & 3. nmean1 differ for each

    consensusToPlot1(groups{Z(i,1)},:) = (consensusToPlot1(groups{Z(i,1)},:))*stdRatio1 + newMean  -mean1*stdRatio1;
    consensusToPlot1(groups{Z(i,2)},:) = (consensusToPlot1(groups{Z(i,2)},:))*stdRatio2 + newMean -mean2*stdRatio2;

    % element
%     consensusToPlot1(groups{Z(i,1)},mask) = (consensusToPlot1(groups{Z(i,1)},mask)-nmean1)*stdRatio1+mean1*meanRatio1;
%     consensusToPlot1(groups{Z(i,2)},mask) = (consensusToPlot1(groups{Z(i,2)},mask)-nmean2)*stdRatio2+mean2*meanRatio2;

%     consensusToPlot1(groups{Z(i,1)},~mask) = (consensusToPlot1(groups{Z(i,1)},~mask)-nmeano1)*stdRatio1+meano1*meanRatio1;
%     consensusToPlot1(groups{Z(i,2)},~mask) = (consensusToPlot1(groups{Z(i,2)},~mask)-nmeano2)*stdRatio2+meano2*meanRatio2;

%     if i > 24
%     figure,imagesc(consensusToPlot1)
%     end
% 
%    
    newVec = nanmean([consensusToPlot1(groups{Z(i,1)},:) ;consensusToPlot1(groups{Z(i,2)},:) ]);
    overlapPx = ~isnan(newVec).*~isnan(outBarsCorrectLengthRescale(1,:));
    newScore(i) = zscore(newVec(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);
        
    oldBarsVec = nanmean([consensusOldMat(groups{Z(i,1)},:) ;consensusOldMat(groups{Z(i,2)},:) ]);

    oldScore(i) = zscore(oldBarsVec(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);

    % old score z-scored


  % old score z-scored
    groups2 = [consensusOldMat(groups{Z(i,1)},:) ;consensusOldMat(groups{Z(i,2)},:) ];
    groups2 = (groups2-nanmean(groups2,2))./nanstd(groups2,1,2);
    oldBarsVec2 = nanmean(groups2);

    oldScore2(i) = zscore(oldBarsVec2(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);


end
%%
figure,plot(newScore);hold on;plot(oldScore);plot(oldScore2);legend({'newScore','oldScore','simple zscore'},'location','southoutside')


figure,plot(consensusToPlot1(curElts(i,1),:));hold on
plot(consensusToPlot1(curElts(i,2),:))
%% test:


% T = cluster(Z,'cutoff',3,'Depth',4);


% all scores 



ix=1
% Start with 1:
bar1 = curSinkVec(ix);
bar2  = curSourceVec(ix);
% figure,imagesc(consensusToPlot1([bar1 bar2],:))
figure,plot(consensusToPlot1([bar1],:));hold on
plot(consensusToPlot1([bar2],:))


b1 = consensusToPlot1(bar1,:);
b2 = consensusToPlot1(bar2,:);


mask =logical( ~isnan(consensusToPlot1(bar1,:)).*~isnan(consensusToPlot1(bar2,:)));

mean1 = mean(b1(mask));
std1 = std(b1(mask),1);

mean2 = mean(b2(mask));
std2 = std(b2(mask),1);

m = nanmean(b1);
s = nanstd(b1,1);
% b1 = (b1-m)./s;

% b1 = 
b2 = (b2-mean2)./std2*std1+mean1;

correctAvg  =[b1;b2];
figure,imagesc([b1;b2])

figure,plot(b1);
hold on
plot(b2)

newVec = nanmean([b1;b2]);
overlapPx = ~isnan(newVec).*~isnan(outBarsCorrectLengthRescale(1,:));
newScore = zscore(newVec(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);
        
oldBarsVec = nanmean([consensusToPlot1(bar1,:);consensusToPlot1(bar2,:)]);

oldScore2 = zscore(oldBarsVec(logical(overlapPx)),1)*zscore(outBarsCorrectLengthRescale(1,logical(overlapPx)),1)'/sum(overlapPx);


%     ii=1;
%     m = mean(barsC{id(ii)}.rawBarcode(barsC{id(ii)}.rawBitmask));
%     s = std(barsC{id(ii)}.rawBarcode(barsC{id(ii)}.rawBitmask),1);
%     barZscored{1}.rawBarcode = (barsC{id(ii)}.rawBarcode-m)./s;
%     barZscored{1}.rawBitmask = barsC{id(ii)}.rawBitmask;
% 
%     barZscored{1}.rawBarcode((barZscored{1}.rawBitmask~=1)) = nan;
%     for ii = 2:length(bars)
% %     dataTable(id(1),1:2)
% 
%         matchStart = dataTable(id(ii),1);
%         matchStop = min(dataTable(id(ii-1),2),dataTable(id(ii),2));
%     
%         b1 =  barZscored{ii-1}.rawBarcode;
%         b1(~barZscored{ii-1}.rawBitmask) = nan;
%     
%         b2 = barsC{id(ii)}.rawBarcode;
%         b2(~barsC{id(ii)}.rawBitmask) = nan;
%     
%         valsb2 = b2(1:matchStop-matchStart+1);
%         valsb1 = b1(matchStart-dataTable(id(ii-1),1):matchStart-dataTable(id(ii-1),1)+length(valsb2)-1);
%     
%         mask = ~isnan(valsb2.*valsb1);
%     
%         mean1 = mean(valsb1(mask));
%         std1 = std(valsb1(mask),1);
%     
%         mean2 = mean(valsb2(mask));
%         std2 = std(valsb2(mask),1); 
% 
%         barZscored{ii}.rawBarcode = (b2-mean2)./std2*std1+mean1;
%         barZscored{ii}.rawBitmask = barsC{id(ii)}.rawBitmask;
%     end




%% Generate data. Need to be able to compare to theory (something like figure S5)

% Schematics of de-novo assembly problem using DNA barcodes.
rng('default')

snr = 2; 
import Core.load_synth_data;
[bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(30,2.72,850,100,snr,0.05,1,2000,150);

sets.minOverlap = 300;
timestamp ='';
idxRun = 1;
oS = synthStr2{idxRun};%resRun{idxRun}.oS;
barcodeGen = bG{idxRun}';
% % 
% % % both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver,pvalCombined,  sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,nan, 0.08,0.01); %localCCStruct

sortedVals = sortedVals(1:length(localScore));
sortedVals(100:end) = nan; % keep less scores to show two islands
% sortedVals(2000:end) = nan; % keep less scores to show two islands
sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.02;
sets.pxDifSetting = 20; %20

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] = ...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, []);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% % % 
% % % for idx = 1:3;
% % %   bars = barcodeGen(barIslands{idx}); 
% % %     barx = barIslands{idx};
% % %    
% % %     print(['FIGS/Fig5S_', num2str(idx), '.eps'],'-depsc','-r300');
% % % 
% % % end
