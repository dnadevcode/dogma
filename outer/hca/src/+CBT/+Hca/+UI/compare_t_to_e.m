function [lm,cache] = compare_t_to_e( lm,ts, cache )
    if nargin < 3
        cache = containers.Map();
    end

    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    flmbs6 = FancyListMgrBtnSet();
    flmbs6.NUM_BUTTON_COLS = 1;
    flmbs6.add_button(compare_theory(ts));

    lm.add_button_sets(flmbs6);
    
  function [btnAddKymos] =compare_theory(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            'COMPARE', ...
            @(~, ~, lm) on_compare_theory(lm, ts));
        
        function [lm] = on_compare_theory(lm, ts)
            
            [selected, ~] = get_selected_list_items(lm);
            
            % change this in the future to support many chromosomes at the
            % same time and do the analysis at the same time
            hcaSessionStruct = cache('hcaSessionStruct');
            hcaSessionStruct.theoryGen.theoryBarcodes = selected(:,2);
            hcaSessionStruct.theoryGen.theoryNames = selected(:,1);
            hcaSessionStruct.theoryGen.bitmask = cellfun(@(x) ones(1,length(x)), hcaSessionStruct.theoryGen.theoryBarcodes,'UniformOutput',false);
                    
            sets = cache('sets');
            sets.barcodeGenSettings = hcaSessionStruct.theoryGen.sets;
             
           % sets.changeBpNmRatio = 1;
            if sets.skipChangeBpNmRatio ~= 1
                titleText = 'Selection of nm/bp ratio (ver.11/12/17)';
                
                import CBT.Hca.UI.get_nmbp_settings;
                [ newNmBp] = get_nmbp_settings(hcaSessionStruct.theoryGen.sets.meanBpExt_nm,titleText);
                
                if newNmBp < hcaSessionStruct.theoryGen.sets.meanBpExt_nm
                    % explain here in a bit more detail how we convert from one
                    % nm/bp ratio to another nm/bp ratio.
                    import CBT.Hca.Core.Analysis.convert_nm_ratio;
                    hcaSessionStruct = convert_nm_ratio(newNmBp,hcaSessionStruct,sets );
                    
                    % make sure that the new mean bp ext is saved for later p-value
                    % calculations
                    sets.barcodeGenSettings.meanBpExt_nm = newNmBp;
                end
            end
                                
             % stretch factor set up
            sets.barcodeConsensusSettings.skipStretch = 1;
            sets.barcodeConsensusSettings.stretchFactors = [1];
            if sets.barcodeConsensusSettings.skipStretch~=0
                import CBT.ExpComparison.Import.prompt_stretch_factors;
                sets.barcodeConsensusSettings.stretchFactors = prompt_stretch_factors();
            end
    
            import CBT.Hca.UI.compare_theory_to_exp;
            hcaSessionStruct = compare_theory_to_exp(hcaSessionStruct, sets);
           

            import CBT.Hca.UI.combine_chromosome_results;
            hcaSessionStruct = combine_chromosome_results(hcaSessionStruct,sets);

            cache('hcaSessionStruct') = hcaSessionStruct;
            cache('sets') = sets;

            import CBT.Hca.UI.launch_export_ui;
            cache = launch_export_ui(ts, cache);
   
            % display results
            import CBT.Hca.UI.get_display_results;
            get_display_results(hcaSessionStruct,sets)
         
            % additional results
            import CBT.Hca.UI.display_additional_results_ui;
            cache = display_additional_results_ui(ts, cache);

            
            import CBT.Hca.UI.display_additional_results_theory_ui;
            cache = display_additional_results_theory_ui(ts, cache);

           % 
        end
    end
end

