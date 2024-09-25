function [cache] = edit_kymographs( lm, ts, cache )
    if nargin < 3
        cache = containers.Map();
    end

    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    flmbs2 = FancyListMgrBtnSet();
    flmbs2.NUM_BUTTON_COLS = 1;
    flmbs2.add_button(edit_kymographs(ts));

    lm.add_button_sets(flmbs2);
    
    function [btnGenerateConsensus] = edit_kymographs(ts)
        % edit_kymographs sets up the structure
        function on_edit_kymographs(lm, ts)
            
            % selected items
            [selectedItems, ~] = get_selected_list_items(lm);
            
            % put the selected items in a structure
            kymoStructs = cell(1,size(selectedItems,1));
            kymoNames =  cell(1,size(selectedItems,1));
            for it=1:size(selectedItems,1)   
            	kymoStructs{it} = selectedItems{it,2};
                kymoNames{it} =  selectedItems{it,1};
            end
            
            % load default settings
            import CBT.Hca.Import.set_default_settings;
            sets = set_default_settings();

            %sets.promptForTimeFr = 1;
            
            % choose timeframes nr for unfiltered kymographs
            if sets.promptForTimeFr ~= 0
                import CBT.Hca.UI.get_hca_settings;
                [ sets.timeFramesNr] = get_hca_settings(sets.timeFramesNr); 
            end
            
            % define session structure
            hcaSessionStruct = struct();

            % put kymos in the session structure
            import CBT.Hca.UI.edit_kymographs_fun;
            hcaSessionStruct = edit_kymographs_fun(hcaSessionStruct,kymoStructs,kymoNames);
            
            % put in the cache for future reference
            cache('hcaSessionStruct') = hcaSessionStruct;     
            cache('sets') = sets;     

        end

        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnGenerateConsensus = FancyListMgrBtn(...
            'Load settings and add kymographs to struture', ...
            @(~, ~, lm) on_edit_kymographs(lm, ts));
    end


end



