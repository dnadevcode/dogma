function [Gtemp,GtempOrig] = create_graph_from_barsetgen(barsetGen, barcodeGen,tsv,plotFig)

    if nargin <4 || isempty(plotFig)
        plotFig = 1;
    end

    G = digraph;
    G = addnode(G,length(barcodeGen));
    G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
    Gundirected = graph;
    Gundirected = addnode(Gundirected,length(barcodeGen));
    Gundirected.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
    
    
    for i=1:length(barsetGen.barcodeIslands)
        for j=1:size(barsetGen.barcodeIslands{i}.edges,1)
            G = addedge(G, barsetGen.barcodeIslands{i}.edges(j,2),  barsetGen.barcodeIslands{i}.edges(j,3), 1);
        end
    end
        %
    Gtemp = G;
    
    listN = 1:length(barcodeGen);
    endNodes = Gtemp.Edges.EndNodes(:);
    curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';
    listN(curNodes) = [];

    GtempOrig = Gtemp;

    Gtemp = rmnode(Gtemp,listN);

    if plotFig
        
        if nargin < 3 || isempty(tsv)
            figure;
        else
            h=tiledlayout(tsv,1,1);
            nexttile(h);
        end
        
        % f = figure
        plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1);
        title('Overlap graph for all barcode islands')

    end

end

