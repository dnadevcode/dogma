function [barcodeConsensusSettings] = get_barcode_params(barcodeConsensusSettings,skipDefaultConsensusSettings)
	
    % check the condition
    if ~skipDefaultConsensusSettings
        % get_barcode_params
        import CBT.Consensus.Import.get_raw_px_edge_length;
        [  barcodeConsensusSettings.prestretchUntrustedEdgeLenUnrounded_pixels,barcodeConsensusSettings.psfSigmaWidth_nm,barcodeConsensusSettings.deltaCut, barcodeConsensusSettings.prestretchPixelWidth_nm] = get_raw_px_edge_length(...
                barcodeConsensusSettings.psfSigmaWidth_nm, ...
                barcodeConsensusSettings.deltaCut, ...
                barcodeConsensusSettings.prestretchPixelWidth_nm,0);            
    end
end