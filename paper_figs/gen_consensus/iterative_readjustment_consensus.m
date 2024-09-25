function [matRep,consensusNew,barCon,stats,cGen] = iterative_readjustment_consensus(aaBarcodeIslandBlockRep,barcodeGen,barIslands, cGenAll,ix,wminC,NN,alphaNu1)

%   iterative_readjustment_consensus
%   

cInitial = cGenAll{ix};

ccthresh = 0;
% NN = 10;
sF = 0.95:0.01:1.05; % max allowed streth/compression
% Compare using pcc
sets.comparisonMethod = 'mass_pcc';
sets.w = nan;


consensusNew = cell(1,NN+1);

consensus = mean(aaBarcodeIslandBlockRep,'omitnan');

consensusNew{1} = consensus;
barCon =[];
barCon{1}.rawBarcode = consensusNew{1};
barCon{1}.rawBitmask = ~isnan(consensusNew{1});

matRep = cell(1,NN);
cGen = cell(1,NN);
matRep{1} = aaBarcodeIslandBlockRep;
for jj = 1:NN
%     jj
    
    % with re-scaling
    barC = barcodeGen(cInitial.idx);
    initialStretch = cellfun(@(x) x.bestBarStretch,cInitial.comparisonStruct);
    
    thr = [];
    thr(1).rawBarcode = consensusNew{jj}(~isnan(consensusNew{jj}));
    thr(1).rawBitmask = ~isnan(thr(1).rawBarcode) ;
    thr(1).isLinearTF = 1;
    
    
    import Core.rescale_barcode_data;
    [barCRescaled] = rescale_barcode_data(barC,sF,initialStretch);
    

    import CBT.Hca.Core.Comparison.hca_compare_distance;
    [rezMaxMP] = hca_compare_distance(barCRescaled, thr, sets );
        
    % create overlap table
    %     info: {'1)PCC, 2) Pos, 3) Or, 4) SecondPos 5) Len '}
    vals = double(zeros(size(rezMaxMP{1}{1},1),4));
    score = zeros(size(rezMaxMP{1}{1},1),1);
    for i=1:length(barCRescaled)
        [maxV,maxPos] = max(rezMaxMP{1}{1}(i,:));
        posStart = double(rezMaxMP{1}{2}(i,maxPos));
        posStop = double(rezMaxMP{1}{2}(i,maxPos)+length(barCRescaled{i}.rescaled{maxPos}.rawBarcode)-1);
        posOr = double(rezMaxMP{1}{3}(i,maxPos));
        curSF =  initialStretch(i)*sF(maxPos);
        score(i) = maxV;
    
    %     if (score(i) < ccthresh) && jj>1 || curSF<sF(1) || curSF > sF(end)
    %         vals(i,:) = valsprev(i,:);% [cInitial.comparisonStruct{i}.pos cInitial.comparisonStruct{i}.pos+cInitial.comparisonStruct{i}.lengthMatch-1 cInitial.comparisonStruct{i}.or cInitial.comparisonStruct{i}.bestBarStretch]; % all these positions are w.r.t. the consensus
    %     else
            vals(i,:) = [posStart posStop posOr curSF]; % all these positions are w.r.t. the consensus
    end

    % keepScores = score>ccthresh;
    keepScores = 1:length(barCRescaled); % keep all
    vals = vals(keepScores,:);
    ids = barIslands{ix}(keepScores);
    
    import Plot.islandsPosStruct;
    cGenNew.comparisonStruct = islandsPosStruct({vals},{ids});
    cGenNew.idx = cInitial.idx;%(score>ccthresh);
    cGenNew.vals = vals;
    try
    import Core.barcode_island_consensus; %barcodeGen,cGenAll, ix, wmin, thresh,alphaNu
    [multiDimBarNew, Z , allscores,clusterBarcodes] = barcode_island_consensus(barcodeGen,{cGenNew}, 1, wminC,[],alphaNu1);
    catch
    % keep previous iteration result.
    end
    length(clusterBarcodes) == length( cInitial.idx);
    consensusNew{jj+1} = mean(multiDimBarNew,'omitnan');
    barCon{jj+1}.rawBarcode = consensusNew{jj+1};
    barCon{jj+1}.rawBitmask = ~isnan(consensusNew{jj+1});
    
    matRep{jj+1} = multiDimBarNew;
%     valsprev = vals;
    cInitial = cGenNew;
    stats.score{jj} = score;
    cGen{jj} = cGenNew;
    cGen{jj}.idx = cGen{jj}.idx(clusterBarcodes);
    cGen{jj}.comparisonStruct = cGen{jj}.comparisonStruct(clusterBarcodes);
    cGen{jj}.vals = cGen{jj}.vals(clusterBarcodes,:);

end


end

