function [] = precompute_pvalue_files(fullPath, params )
    % precompute_pvalue_files
    
    % Pregenerates p-value database file
    %     Args:
    %         fullPath, params 
    % 
    %     Returns:
    % 
    %     Example:
    %
        
    % compute the psf
    psf = params.psfSigmaWidth_nm/params.pixelWidth_nm;
    
    % compute the long random barcode
    rand2 = normrnd(0, 1, 1, params.len2);
    % convolve with a Gaussian
    import CBT.Hca.Core.Pvalue.convolve_bar;
    rand2 = convolve_bar(rand2, psf, length(rand2));
    
    % the short lengths that we need to compute for
    import CBT.Hca.Core.Pvalue.compute_random_max_cc;
	import CBT.Hca.Import.load_pval_struct;

    % lengths to compute
    lens =params.lenMin:params.lenMax;
    
    
    [ vals, data ] = load_pval_struct(fullPath);
    % 
    for i =1:length(lens)
        % we load the data every time, maybe change this with simpler
        % approach where we save the data only in the end
        
        lenCur = lens(i);
        disp(strcat(['Computing p-value for barcodes of length ' num2str(lenCur) ', already done ' num2str(i-1) ' out of ' num2str(length(lens))]));
        [dataNew] = compute_random_max_cc(lenCur,rand2, psf, params.numRnd);%%    

        % checks if this length was already considered before         
        [a,b] = ismember(lenCur, vals);
        if a~=0 
            data{b} =[data{b}; dataNew'];
        else
            vals(end+1) = lenCur;
            data{end+1} = dataNew';
        end
    end
    
    % exports structure
    import CBT.Hca.Export.export_pval_struct;
    export_pval_struct( fullPath,vals,data );
        
end

