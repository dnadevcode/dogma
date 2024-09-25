function [lm] = launch_theory_ui(ts)
    % launch_theory_ui -
    %   adds a tab with list management UI/functionality 
    %    
	tabTitle = 'Theory';
    [hTabTheoryImport] = ts.create_tab(tabTitle);
    hPanelTheoryImport = uipanel(hTabTheoryImport);
    ts.select_tab(hTabTheoryImport);

	import Fancy.UI.FancyList.FancyListMgr;
    lm = FancyListMgr();
    lm.set_ui_parent(hPanelTheoryImport);
    lm.make_ui_items_listbox();
    
    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    
    flmbs1 = FancyListMgrBtnSet();
    flmbs1.NUM_BUTTON_COLS = 2;
    flmbs1.add_button(make_add_sequences(ts));
    flmbs1.add_button(make_remove_consensus_btn());   
    
    flmbs2 = FancyListMgrBtnSet();
    flmbs2.NUM_BUTTON_COLS = 2;
	flmbs2.add_button(FancyListMgr.make_select_all_button_template());
    flmbs2.add_button(FancyListMgr.make_deselect_all_button_template());

    function [btnAddKymos] =make_add_sequences(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            'Add sequence(s)', ...
            @(~, ~, lm) on_add_sequences_directly(lm, ts));
            
        function [] = on_add_sequences_directly(lm, ts)
            [barcodeFilenames, barcodeFilenamesDirpath] = uigetfile(...
                {'*.fasta,*.fa'}, ...
                'Select sequence(s) to import', ...
                pwd, ...
                'MultiSelect','on');
                    
            if ~iscell(barcodeFilenames)
                barcodeFilenames = {barcodeFilenames};
            end
            numFiles = size(barcodeFilenames,2);

            sequencePath = cell(1,numFiles);
            nameSequence = cell(1,numFiles);
            for i=1:numFiles
                sequencePath{i} = strcat(barcodeFilenamesDirpath,barcodeFilenames{i});
                nameSequence{i} = barcodeFilenames{i};
            end

            lm.add_list_items(nameSequence, sequencePath);
            end
    end

   
    lm.add_button_sets(flmbs1,flmbs2);

    function [btnRemoveConsensus] = make_remove_consensus_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnRemoveConsensus = FancyListMgrBtn(...
            'Remove selected sequence(s)', ...
            @(~, ~, lm) on_remove_selected_consensus(lm));
        function [] = on_remove_selected_consensus(lm)
            lm.remove_selected_items();
        end
    end

end