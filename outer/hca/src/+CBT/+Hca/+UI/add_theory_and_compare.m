function [cache] = add_theory_and_compare( lm, ts, cache )
    if nargin < 3
        cache = containers.Map();
    end

    % creates a button set
    import Fancy.UI.Templates.create_button_set_ts;
    create_button_set_ts(lm,ts, @make_theory);
    
    function [btnOut] = make_theory(ts)   
        function on_make_theory(lm, ts)
            
            tabTitle = 'Theory';
            [hTabTheoryImport] = ts.create_tab(tabTitle);
            hPanelTheoryImport = uipanel(hTabTheoryImport);
            ts.select_tab(hTabTheoryImport);

            import Fancy.UI.FancyList.FancyListMgr;
            lm = FancyListMgr();
            lm.set_ui_parent(hPanelTheoryImport);
            lm.make_ui_items_listbox();
            
            import CBT.Hca.UI.load_theory_ui;
            [lm,cache] = load_theory_ui(lm,ts,cache);
        
%             import CBT.Hca.UI.compute_theory_ui;
%             [lm,cache] = compute_theory_ui(lm,ts,cache);

            import CBT.Hca.UI.compare_t_to_e;
            [lm,cache] = compare_t_to_e(lm,ts,cache);
            
        end

        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnOut = FancyListMgrBtn(...
            'Add theory and compare', ...
            @(~, ~, lm) on_make_theory(lm, ts));
    end
end

