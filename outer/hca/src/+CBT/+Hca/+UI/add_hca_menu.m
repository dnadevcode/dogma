function [] = add_hca_menu(hMenuParent, tsHCA)
    % add_hca_menu adds menu functions to the table
    %
    % :param hMenuParent: gui window
    % :param tsHCA: gui window

    % :returns: 
        
    % 
    hMenuETE = uimenu('Parent', hMenuParent,'Label', 'HCA');
    
    % Fragment to theory comparison
    import CBT.Hca.UI.run_fragment_vs_theory_similarity_analysis;
    uimenu(hMenuETE,'Label', 'Analyze Fragment vs. Theory Similarity', 'Callback', @(~, ~) run_fragment_vs_theory_similarity_analysis(tsHCA));
   
    % loading results that there calculatated by run_fragment_vs_theory_similarity_analysis
    import CBT.Hca.Import.load_hca_results;
    uimenu(hMenuETE,'Label', 'Load HCC Results', 'Callback', @(~, ~) load_hca_results(tsHCA));

   
end