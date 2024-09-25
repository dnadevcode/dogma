function [barcodeGenData] = gen_barcode_data(alignedKymo,leftEdgeIdxs, rightEdgeIdxs, simple)
    % gen_barcode_data
    %
    %     Args:
    %         alignedKymo: aligned kymograph
    %         leftEdgeIdxs: left edge indices on the aligned kymograph
    %         rightEdgeIdxs: right edge indices on the aligned kymograph
    % 
    %     Returns:
    %         barcodeGenData: Return structure
    % 

    if nargin < 4
        simple = 0;
    end
    
    if simple == 1
        barcodeGenData.rawBarcode = nanmean(alignedKymo,1);
        barcodeGenData.bgMeanApprox = nan;
        barcodeGenData.bgStdApprox = nan;
        barcodeGenData.lE = 1;
        barcodeGenData.rE = length(barcodeGenData.rawBarcode);
        barcodeGenData.rawBg = nan;
    else
        barcodeGenData.leftEdgeIdxs = leftEdgeIdxs;
        barcodeGenData.rightEdgeIdxs = rightEdgeIdxs;
        
        leftEdgeIdx = round(nanmean(leftEdgeIdxs));
        rightEdgeIdx = round(nanmean(rightEdgeIdxs));

        if rightEdgeIdx <= leftEdgeIdx
            disp('ah');
        end
        
        % Determine indices for rotated barcode with background cropped out
        barcodeIdxs = leftEdgeIdx:rightEdgeIdx;
        rawBarcode = nanmean(alignedKymo, 1);

        % define background indices
        bgIndices = true(1, size(alignedKymo, 2));
        bgIndices(barcodeIdxs) = 0;

        nonBarcodeVals = alignedKymo(:, bgIndices);
        nonBarcodeVals = nonBarcodeVals(:);
        
        if isempty(nonBarcodeVals)
            barcodeGenData.bgMeanApprox =  nan;
            barcodeGenData.bgStdApprox = nan;
        else
            nonBarcodeDistFit = fitdist(nonBarcodeVals, 'Normal');
            barcodeGenData.bgMeanApprox =  nonBarcodeDistFit.mu;
            barcodeGenData.bgStdApprox = nonBarcodeDistFit.sigma;
        end


        rawBg = nanmean(rawBarcode(bgIndices));

        barcodeGenData.rawBarcode = rawBarcode(barcodeIdxs);

        barcodeGenData.lE = leftEdgeIdx;
        barcodeGenData.rE = rightEdgeIdx;
        barcodeGenData.rawBg = rawBg;
    end
end
