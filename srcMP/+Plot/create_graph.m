

function G=create_graph(barcodeGen,PCC_OVERLAP,sortedIds,NN)
    % first create barcode graph from all:
    G = digraph;
    G = addnode(G,length(barcodeGen));
    G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
    nodeOrientations = ones(1,length(barcodeGen));
    for ii=1:NN
        [curSink,curSource] = ind2sub(size(PCC_OVERLAP),sortedIds(ii));
        G = addedge(G, curSource, curSink, 1); 
    end
    Gtemp = G;
    % %         
    listN = 1:length(barcodeGen);
    endNodes = Gtemp.Edges.EndNodes(:);
    curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';
    listN(curNodes) = [];
    % %     
    Gtemp = rmnode(Gtemp,listN);
    % %     weak_bins = conncomp(Gtemp,'Type','weak');
    % %         numComp(idxPair) = max(weak_bins);
    % %         numNodes(idxPair) = length(curNodes);
    % %     
    % %     
    figure;
    plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1)
% %     

end