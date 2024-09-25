function [lm,cache] = compute_theory_ui(lm,ts,cache)
  % compute theory UI
    if nargin < 3   % but no reason this should be the case
        cache = containers.Map();
    end

    import Fancy.UI.FancyList.FancyListMgrBtnSet; 
    flmbs5 = FancyListMgrBtnSet();

    flmbs5.NUM_BUTTON_COLS = 1;

    flmbs5.add_button(select_theory(ts));

     function [btnAddKymos] =select_theory(ts)
            import Fancy.UI.FancyList.FancyListMgrBtn;
            btnAddKymos = FancyListMgrBtn(...
                'Compute theory barcodes for sequences', ...
                @(~, ~, lm) on_select_theory(lm, ts));


            function [] = on_select_theory(lm, ts)
                [selectedItems, ~] = get_selected_list_items(lm);

                % load settings
                sets = cache('sets');
                if ~sets.skipBarcodeGenSettings
                    import CBT.Hca.UI.comparison_settings;
                    sets.barcodeGenSettings = comparison_settings(); %
                end

                import CBT.Hca.Core.Theory.compute_free_concentrations;
                sets.barcodeGenSettings = compute_free_concentrations(sets.barcodeGenSettings);

                for i=1:size(selectedItems,1)
                    FASTAData = fastaread(selectedItems{i,2});
                    seq = FASTAData(1).Sequence;
                    name = FASTAData(1).Header;
                    disp(strcat(['loaded theory sequence ' name] ));

                    import CBT.Hca.Core.compute_hca_theory_barcode;
                    [theorySeq,bitmask] = compute_hca_theory_barcode(seq,sets.barcodeGenSettings);
                    theoryGen.theoryBarcodes{i} = theorySeq;
                    theoryGen.theoryNames{i} = name;
                    theoryGen.bitmask{i} = bitmask;
                    theoryGen.bpNm{i} = sets.barcodeGenSettings.meanBpExt_nm/sets.barcodeGenSettings.psfSigmaWidth_nm;

                 end
                    theoryGen.sets = sets.barcodeGenSettings;

                    hcaSessionStruct = cache('hcaSessionStruct');
                    hcaSessionStruct.theoryGen = theoryGen;
                    cache('hcaSessionStruct') = hcaSessionStruct;
                    cache('sets') = sets;
                end
     end
     lm.add_button_sets(flmbs5);
end
