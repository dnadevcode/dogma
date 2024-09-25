function [kymoStructs,barGen] = load_kymo_data(sets)

    sets.kymoImportMethod = 'tifs';
    switch sets.kymoImportMethod
        case 'tifs'
            import CBT.Hca.Import.add_kymographs_fun;
            [kymoStructs] = add_kymographs_fun(sets);

            import CBT.Hca.Core.edit_kymographs_fun;
            kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr);



            % align kymos
            import CBT.Hca.Core.align_kymos;
            [kymoStructs] = align_kymos(sets,kymoStructs);
            
            % generate barcodes
            import CBT.Hca.Core.gen_barcodes;
            barGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);
        otherwise
    end

end

