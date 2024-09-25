function [lm,cache] = load_theory_ui(lm,ts,cache)
     if nargin < 3
        cache = containers.Map();
     end

    import Fancy.UI.FancyList.FancyListMgrBtnSet;
	import Fancy.UI.FancyList.FancyListMgr;
    flmbs2 = FancyListMgrBtnSet();
    flmbs2.NUM_BUTTON_COLS = 2;
    flmbs2.add_button(FancyListMgr.make_select_all_button_template());
    flmbs2.add_button(FancyListMgr.make_deselect_all_button_template());
    lm.add_button_sets(flmbs2);

        
    
    flmbs3 = FancyListMgrBtnSet();

    flmbs3.NUM_BUTTON_COLS = 2;

    flmbs3.add_button(make_add_sequences(ts));
    flmbs3.add_button(make_remove_consensus_btn());   

    
    

    function [btnAddKymos] =make_add_sequences(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            'Add precalculated theory (HCA_theory)', ...
            @(~, ~, lm) on_add_sequences_directly(lm, ts));
        
        
        function [] = on_add_sequences_directly(lm, ts)

            [barcodeFilenames, barcodeFilenamesDirpath] = uigetfile(...
                {'*.mat;'}, ...
                'Select precalculated theory barcodes to import', ...
                pwd, ...
                'MultiSelect','on');
            
            filePath =strcat(barcodeFilenamesDirpath,barcodeFilenames);
            
            hcaSessionStruct = cache('hcaSessionStruct');

            % this is a bit messy still..
            import CBT.Hca.UI.Helper.load_theory;
            hcaSessionStruct = load_theory( filePath,hcaSessionStruct );
            
            cache('hcaSessionStruct') = hcaSessionStruct;
             
            % add only name instead of the whole file for speed!
            lm.add_list_items(hcaSessionStruct.theoryGen.theoryNames, hcaSessionStruct.theoryGen.theoryBarcodes);
            end
    end

    function [btnRemoveConsensus] = make_remove_consensus_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnRemoveConsensus = FancyListMgrBtn(...
            'Remove selected sequence(s)', ...
            @(~, ~, lm) on_remove_selected_consensus(lm));
        function [] = on_remove_selected_consensus(lm)
            lm.remove_selected_items();
        end
    end

    lm.add_button_sets(flmbs2,flmbs3);


end