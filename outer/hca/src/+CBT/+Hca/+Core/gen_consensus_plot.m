function [  ] = gen_consensus_plot( consensusStructs,consensusStruct,barcodeGen, sets )
    %gen_consensus_plot
    
%	load('consensus_struct.mat');

          % loads figure window
    hFig = figure(...
        'Name', 'CB HCA tool', ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 1 1], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none' ...
    );
    hMenuParent = hFig;
    hPanel = uipanel('Parent', hFig);
    import Fancy.UI.FancyTabs.TabbedScreen;
    ts = TabbedScreen(hPanel);
    
    import CBT.Hca.Export.plot_consensus_concentrically;
%     consensusStructs,barcodeGen, idx, ts
    plot_consensus_concentrically(consensusStructs,barcodeGen,1,ts)
    
    %%
%         import CBT.Consensus.UI.plot_clusters_linearly;
%     plot_clusters_linearly(ts, consensusStruct);
    
    cs = clusterConsensusData.details.consensusStruct;
           % loads figure window
    hFig = figure(...
        'Name', 'CB HCA tool', ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 1 1], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none' ...
    );
    hMenuParent = hFig;
    hPanel = uipanel('Parent', hFig);
    import Fancy.UI.FancyTabs.TabbedScreen;
    ts = TabbedScreen(hPanel);
    
    import CBT.Consensus.UI.generate_consensus_ui;
    generate_consensus_ui(ts, cs);
    
    %%
    import CBT.Consensus.UI.plot_clusters_concentrically;
    hTabConcentric = plot_clusters_concentrically(ts, clusterKeys, clusterConsensusDataStructs);

    
    
 
%     import CBT.Consensus.UI.launch_export_ui;
%     launch_export_ui(ts, clusterKeys, clusterConsensusDataStructs)

    
    
    
    
    import CBT.Consensus.UI.generate_consensus_ui;
    generate_consensus_ui(ts, consensusStruct);
    
    import CBT.Consensus.UI.plot_clusters_linearly;
    plot_clusters_linearly(ts, consensusStruct);

    hTabPairwiseConsensusHistory = ts.create_tab('Pairwise Consensus History');
    ts.select_tab(hTabPairwiseConsensusHistory);
    hPanelPairwiseConsensusHistory = uipanel(hTabPairwiseConsensusHistory);
    import CBT.Consensus.UI.plot_pairwise_consensus_history;
    plot_pairwise_consensus_history(consensusStruct, hPanelPairwiseConsensusHistory);

    hTabConsensusDendros = ts.create_tab('Consensus Dendrograms');
    ts.select_tab(hTabConsensusDendros);
    hPanelConsensusDendros = uipanel(hTabConsensusDendros);
    import CBT.Consensus.UI.plot_consensus_dendrograms;
    plot_consensus_dendrograms(consensusStruct, hPanelConsensusDendros);

    
end

