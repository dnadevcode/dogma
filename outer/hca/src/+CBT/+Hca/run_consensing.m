function [] = run_consensing(tsCBC)
    tabTitle = 'Kymos for Consensus';
    hTabKymoImport = tsCBC.create_tab(tabTitle);
    tsCBC.select_tab(hTabKymoImport);
    hPanelKymoImport = uipanel(hTabKymoImport);

    import CBT.Consensus.UI.launch_kymo_import_ui;
    lm = launch_kymo_import_ui(hPanelKymoImport, tsCBC);

    import CBT.Consensus.UI.add_consensus_btns_to_kymo_list_mgr;
    add_consensus_btns_to_kymo_list_mgr(lm, tsCBC);
end