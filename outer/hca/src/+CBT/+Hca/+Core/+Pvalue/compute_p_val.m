function [ pvalResults ] = compute_p_val(pvalData, comparisonStruct, barcodeGen, consensusStruct, sets )
    % compute_p_val
    %
    % Compute p-values for comparison results from before
    %     Args:
    %         pvalData, comparisonStruct
    % 
    %     Returns:
    %         pvalResults  
    %     Example:
    %
    % maximum p-values
    
    cMaxVals = cellfun(@(x) x.maxcoef(1), comparisonStruct);
    
    % length of theory barcodes
    barLen = cellfun(@(x) sum(x.rawBitmask), barcodeGen);
    
    %% if there is a consensus
    if sets.genConsensus == 1
       barLen(end+1) = sum(consensusStruct.bitmask);
    end
    
    % stretch factors
    strFac = sets.theory.stretchFactors;
    
    import CBT.Hca.Core.Pvalue.compute_p_val_score;
    [ p, cal ] = compute_p_val_score(cMaxVals, pvalData, barLen,strFac );

    pvalResults.pValueMatrix = p;
    pvalResults.pValueCalculated = cal;              
end

