% we split data into two to generate two consensus

ix=2;
cGenTemp = [];
cGenTemp{1} = cGenAll{ix};
cGenTemp{1}.idx  = cGenTemp{1}.idx(1:end/2);
cGenTemp{1}.comparisonStruct = cGenTemp{1}.comparisonStruct(1:end/2);
cGenTemp{2} = cGenAll{ix};
cGenTemp{2}.idx  = cGenTemp{2}.idx(end/2+1:end);
cGenTemp{2}.comparisonStruct = cGenTemp{2}.comparisonStruct(end/2+1:end);

barIslandsTemp{1} = barIslands{ix}(1:end/2);
barIslandsTemp{2} = barIslands{ix}(end/2+1:end);

wminC = 300;
NN = 5;

import Core.barcode_island_consensus;
aaBarcodeIslandBlockRep = cell(1,length(cGenAll));
clusterBarcodes =  cell(1,length(cGenAll));
matRep =  cell(1,length(cGenAll));
consensusNew =  cell(1,length(cGenAll));
barCon =  cell(1,length(cGenAll));
stats =  cell(1,length(cGenAll));
cIt =  cell(1,length(cGenAll));

for ix = 1:length(cGenTemp)
%     ix
    [aaBarcodeIslandBlockRep{ix}, ~,~,clusterBarcodes{ix},~] = barcode_island_consensus(barcodeGen,cGenTemp, ix, wminC);

    NN = 5;
    % wminC
    [matRep{ix},consensusNew{ix},barCon{ix},stats{ix},cIt{ix}] = iterative_readjustment_consensus(aaBarcodeIslandBlockRep{ix},barcodeGen, barIslandsTemp, cGenTemp,ix,wminC,NN);

end


outConsensus = cellfun(@(x) x{end},matRep,'un',false);

