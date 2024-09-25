function [goodPos,passLocalThresh,passGlobalThresh,sortedVals,sortedIdsAll,foS] = get_pairs_passing_thresh(oS,mpMaxLenBased,minOverlap,minLen,Nvals)
%  find which barcodes pass the thresh
    [sortedVals, sortedIdsAll, pscores,foS,overlaplen] = sorted_scores(oS);
%     posShift = calc_pos_dif(oS,sortedIdsAll,synthStr{idxRun},lenThry, Nvals);
    
    import Nullmodel.filter_scores_local;
    import Nullmodel.filter_scores_global;
    [passLocalThresh] = filter_scores_local(sortedVals,mpMaxLenBased,minOverlap,minLen,Nvals);
    [passGlobalThresh] = filter_scores_global(foS,mpMaxLenBased,overlaplen,minLen,Nvals);
    
    goodPos = passLocalThresh.*passGlobalThresh;
end

