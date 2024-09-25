function [commonLength,aborted] = make_barcode_settings(rawBarcodeLens)

    commonLength = ceil(mean(rawBarcodeLens));

    import CBT.Consensus.Import.confirm_stretching_is_ok;
    [notOK, commonLength] = confirm_stretching_is_ok(commonLength, rawBarcodeLens);
    aborted = notOK;
    if aborted
        fprintf('Skipping consensus generation\n');
        return;
    end
   
    if aborted
        fprintf('Aborting consensus input generation\n');
        return;
    end
    
end