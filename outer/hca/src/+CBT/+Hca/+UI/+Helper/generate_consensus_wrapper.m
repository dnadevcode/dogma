function [ consensusStruct,rawBarcodes,rawBitmasks,displayNames ] = generate_consensus_wrapper(rawBarcodes,rawBitmasks,barcodeConsensusSettings,barcodeGen,displayNames,lens, cache )
         
        if nargin < 7
            cache = containers.Map();
        end

        rawBgs = cellfun(@(xx) xx.bgMeanApprox,barcodeGen,'UniformOutput',false);

        import CBT.Consensus.Core.generate_consensus_for_barcodes;
        [consensusStruct, cache] = generate_consensus_for_barcodes(rawBarcodes, displayNames,-1*ones(length(rawBarcodes),1),barcodeConsensusSettings, cache, rawBgs);
      % add consensus barcode as the last barcode in the structure
        lengths = cellfun(@length,consensusStruct.clusterKeys);
        [~,b] = max(lengths);
        key = consensusStruct.clusterKeys{b};
        consSt = consensusStruct.barcodeStructsMap(key);
        
        displayNames{length(lens)+1} = key ;
        rawBarcodes{length(lens)+1} = consSt.barcode;        
        rawBitmasks{length(lens)+1} = logical(consSt.indexWeights);
     
                
      
end

