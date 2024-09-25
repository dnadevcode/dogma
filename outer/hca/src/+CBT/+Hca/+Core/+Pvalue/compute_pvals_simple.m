function [pval] = compute_pvals_simple(theoryStruct,comparisonStruct,sets)

    import CBT.Hca.Core.Pvalue.compute_pvals_simple;

    % If one computes these directly from data, then we want a random
        % theory and 1000 random barcodes. We then compare compare these (
        % so takes 1000 times) using compare distance. Finally we compute
        % p-values
        sets.circ = 0;
        sets.kernelSigma = theoryStruct{1}.psfSigmaWidth_nm/theoryStruct{1}.pixelWidth_nm; % same as first theory
        sets.length1 = max(cellfun(@(x) x.bestLength,comparisonStruct));
        sets.theory.stretchFactors = 1; % for pval, only for the best
%         tic
         p = @(x,par1,par2) 1-(0.5+0.5*(1-betainc((x).^2,0.5,par1/2-1,'upper'))).^par2 ;

        import CBT.Hca.Import.load_data_wrapper;
        import CBT.Hca.Core.Pvalue.compute_evd_params;
         import CBT.Hca.Core.Comparison.compare_distance;

%         [~,~,barcodeGenRand] = load_data_wrapper('random',sets.output.matDirpath,sets,sets.kymosets.kymofilefold{1} )
%     toc

        pval = ones(1,length(comparisonStruct));
        % now loop through these
        for i=1:length(comparisonStruct)
            sets.length1 = comparisonStruct{i}.bestLength;  
            if sets.length1 >= sets.w
                [~,~,barcodeGenRand] = load_data_wrapper('random',[],sets,[] );
                theoryStructP = {theoryStruct{comparisonStruct{i}.idx}};
                [rezMax,~,~] = compare_distance(barcodeGenRand,theoryStructP, sets, [] );
                %
                coeffs = cellfun(@(x) x.maxcoef,rezMax{1});
                %
                params = compute_evd_params(coeffs,100);
                pval(i) = p(comparisonStruct{i}.maxcoef(1),params(1),params(2));
            else
                pval(i) = 1;
            end
        end
        
%             timestamp = datestr(clock(),'yyyy-mm-dd HH:MM:SS');

        
    fid = fopen(fullfile(sets.output.matDirpath,sets.timestamp,'pval_results.txt'), 'w');
    fprintf(fid, '#Date: %s\n', sets.timestamp);
    fprintf(fid, '%6.5e ', pval);
    % maybe add some of the latter
%     fprintf(fid, '#Name of barcode 1: %s\n', barcodeNameA);
%     fprintf(fid, '#Name of barcode 2: %s\n', barcodeNameB);
%     fprintf(fid, '#Stretched to same length: %s\n', sameLengthTFStr);
%     fprintf(fid, '#Barcode 1 stretched for matching: %.0f%%\n', barcodeStretchPercentageA);
%     fprintf(fid, '#Kbp per pixel: %s\n', num2str(kbpsPerPixel));
%     fprintf(fid, '#P-value: %g\n', pValAB);
%     fprintf(fid, '#Cross-correlation: %g\n', ccAB);
%     fprintf(fid, '#\n');

    % Write first row (containing barcode names)
%     fprintf(fid, '#Intensity 1\t\tIntensity 2\t\t\n');

%     % Write the rest
%     maxPixelIdx = length(longerBarcodeShifted);
%     for pixelIdx = 1:maxPixelIdx
%         fprintf(fid, '%5.4f\t\t\t\t%g\n', barcodeShiftedA(pixelIdx), barcodeShiftedB(pixelIdx));
%     end
    fclose(fid);
        
        % now save in output folder
end

