function [] = run_fragment_vs_theory_similarity_analysis(tsHCC)
    % run_fragment_vs_theory_similarity_analysis
    % Takes kymographs as an input and computes the consensus
    %
    %     Args:
    %         tsHCC (figure): Input tab
    % 
    %     Returns:
    % 
    %     Example:
    %          run_fragment_vs_theory_similarity_analysis(tsHCC)


    % create a main tab for analysis and importing kymographs:
	import Fancy.UI.Templates.create_import_tab;
    lm = create_import_tab(tsHCC,'Fragment vs theory similarity analysis');
  
    % load settings and kymo structure
    import CBT.Hca.UI.load_settings_and_session_structure;
    cache = load_settings_and_session_structure(lm);

    % make consensus structure
    import CBT.Hca.UI.make_consensus_structure;
    cache = make_consensus_structure(lm, cache);

    % export structure with barcodes and consensus
    import CBT.Hca.UI.export_consensus_structure;
    cache = export_consensus_structure(lm, cache);
    
    % add theory and compare
    import CBT.Hca.UI.add_theory_and_compare;
    cache = add_theory_and_compare(lm, tsHCC,cache);

end