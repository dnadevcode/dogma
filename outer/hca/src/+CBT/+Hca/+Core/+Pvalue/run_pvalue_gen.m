function [ pvalResults ] = run_pvalue_gen( comparisonStruct, barcodeGen, consensusStruct, sets )
    % run_pvalue_gen

    import CBT.Hca.Import.load_pval_struct;
    import CBT.Hca.Core.Pvalue.pregenerate_pvalue_db;

    if sets.pvalue.generatemoreccs == 1
         pregenerate_pvalue_db(sets.pvalue);
    end
    
    try 
        fullPath = strcat([sets.pvalue.fold,sets.pvalue.file]);  
        [ pvalData.len1, pvalData.data ] = load_pval_struct(fullPath);
    catch
         disp('No pre-computed p-value database chosen. Running precompute method... ');
         % make want to make cure pvalue folder exists?
         sets.pvalue.file = 'pval.txt';
         pregenerate_pvalue_db(sets.pvalue);
         fullPath = strcat([sets.pvalue.fold,sets.pvalue.file]);  
         [ pvalData.len1, pvalData.data ] = load_pval_struct(fullPath);
    end
            
    for i=1:length(pvalData.data)
        pvalData.data{i} = pvalData.data{i}(1:sets.pvalue.numRnd);
    end
    import CBT.Hca.Core.Pvalue.compute_p_val;
    [ pvalResults ] = compute_p_val(pvalData, comparisonStruct, barcodeGen, consensusStruct, sets);
         

end


