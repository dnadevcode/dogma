function [] = load_hca_results( tsHCA, hcaSessionStruct )
   if nargin < 2
        import CBT.Hca.Import.load_hca_struct;
        hcaSessionStruct = load_hca_struct();
        sets = hcaSessionStruct.sets;
        if isempty(hcaSessionStruct)
            return;
        end
    end
    if isempty(hcaSessionStruct)
        return;
    end



   % assignin('base', 'hcaSessionStruct', hcaSessionStruct);
    import CBT.Hca.UI.get_display_results;

    cache = containers.Map();
    cache('hcaSessionStruct') = hcaSessionStruct;
    cache('sets') = sets;

    get_display_results(hcaSessionStruct,sets)
    fprintf('Finished loading HCA session results \n');

    import CBT.Hca.UI.display_additional_results_ui;
    display_additional_results_ui(tsHCA, cache);
    
    import CBT.Hca.UI.display_additional_results_theory_ui;
    display_additional_results_theory_ui(tsHCA, cache);
    
    import CBT.Hca.UI.launch_export_ui;
    launch_export_ui(tsHCA, cache);
            

end

