function [barConRef,Z , allscores,clusterBarcodes] = ref_based_consensus(theoryStruct,cGenAll,ix,barcodeGen,barIslands,wminC,sFRe)


thrRef = [];
thrRef(1).rawBarcode = theoryStruct.rawBarcode;
thrRef(1).rawBitmask = ~isnan(thrRef(1).rawBarcode) ;
thrRef(1).isLinearTF = 1;

barC = barcodeGen( cGenAll{ix}.idx);

initialStretch  = ones(1,length(barC));

import Core.rescale_barcode_data;
[barCRescaled] = rescale_barcode_data(barC,sFRe,initialStretch);

sets.comparisonMethod = 'mass_pcc';
sets.w = nan;
import CBT.Hca.Core.Comparison.hca_compare_distance;
[rezMaxMP] = hca_compare_distance(barCRescaled, thrRef, sets );
    
% create overlap table
%     info: {'1)PCC, 2) Pos, 3) Or, 4) SecondPos 5) Len '}
valsRefbased = double(zeros(size(rezMaxMP{1}{1},1),4));
score = zeros(size(rezMaxMP{1}{1},1),1);
for i=1:length(barCRescaled)
    [maxV,maxPos] = max(rezMaxMP{1}{1}(i,:));
    posStart = double(rezMaxMP{1}{2}(i,maxPos));
    posStop = double(rezMaxMP{1}{2}(i,maxPos)+length(barCRescaled{i}.rescaled{maxPos}.rawBarcode)-1);
    posOr = double(rezMaxMP{1}{3}(i,maxPos));
    curSF =  initialStretch(i)*sFRe(maxPos);
    score(i) = maxV;

    valsRefbased(i,:) = [posStart posStop posOr curSF]; % all these positions are w.r.t. the consensus

%     stats.posdif(i) = cGenAll{ix}.comparisonStruct{i}.pos-posStart;
%     stats.ordif(i) = isequal(cGenAll{ix}.comparisonStruct{i}.or,posOr);
%     stats.sfDif(i) = sF(maxPos);
end

    % keepScores = score>ccthresh;
    keepScores = 1:length(barCRescaled); % keep all
    valsRefbased = valsRefbased(keepScores,:);
    ids = barIslands{ix}(keepScores);
    
    import Plot.islandsPosStruct;
    cGenNewRef.comparisonStruct = islandsPosStruct({valsRefbased},{ids});
    cGenNewRef.idx = ids;%(score>ccthresh);
    import Core.barcode_island_consensus;
    [multiDimBarNew,Z , allscores,clusterBarcodes]  = barcode_island_consensus(barcodeGen,{cGenNewRef}, 1, wminC);
    
    consensusNenRef = mean(multiDimBarNew,'omitnan');
    barConRef{1}.rawBarcode = consensusNenRef;
    barConRef{1}.rawBitmask = ~isnan(consensusNenRef);


end

