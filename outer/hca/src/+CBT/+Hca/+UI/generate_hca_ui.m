function [] = generate_hca_ui(ts, hcaSessionStruct)
  
    import CBT.ExpComparison.UI.struct2vars; %todo: move this somewhere?
    struct2vars(hcaSessionStruct)

    % Create tab
    hTabHCA = ts.create_tab('Human chromosome analysis results');
    ts.select_tab(hTabHCA);
    hPanelHCA = uipanel('Parent', hTabHCA);
    
end