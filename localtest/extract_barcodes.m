function [bG] = extract_barcodes(kymoStructs)

sets.maxLen=Inf;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 6; % depending on nm/bp
sets.minLen = 300; % dependent on nm/px ratio

%     % sets.minLen
bG = cell(1,length(kymoStructs) );
for idFold = 1:length(kymoStructs)

    import DBM4.gen_barcodes_from_kymo;
    barcodeGen =  gen_barcodes_from_kymo( kymoStructs{idFold}, sets,sets.maxLen);
    
    % bitmask by
    % 1) filter out periodic barcodes
    % 2) filter out sigmoid
%     sigm = @(x) 1./(1+exp(-x));
%     NN = 30;
%     filtSigmoid = sigm(-NN:1:NN);
%     import SignalRegistration.masked_pcc_corr;
    % run unmasked comparison
%     for ix=1:length(barcodeGen)
%         bLong =  bars{ix}.rawBarcode(  bars{ix}.rawBitmask);
%         % regular MASS PCC
%         [dist,numElts] = masked_pcc_corr(filtSigmoid,[bLong zeros(1,length(filtSigmoid))],ones(1,length(filtSigmoid)),[ones(1,length(bLong)) zeros(1,length(filtSigmoid))]);
%     figure,plot(dist')
%     end


%     for i=1:length(barcodeGen)
%         barcodeGen{i}.rawBitmask =  barcodeGen{i}.rawBitmask &...
%             abs(barcodeGen{i}.rawBarcode - nanmean(barcodeGen{i}.rawBarcode(barcodeGen{i}.rawBitmask))) ...
%             <= 2 * nanstd(barcodeGen{i}.rawBarcode(barcodeGen{i}.rawBitmask));
%     end

    bG{idFold} = barcodeGen;
end

end

