function [finalgraph,Ggraphs,false_maps] = create_overlap_graph(sizeOs,sortedIds,NN,barcodeGen,curStd,thrStd,overlapStruct,tsv)
    
    % Overlap graph - visualization of bargrousp
    %   Args:
    %
    %   Returns:
    %       finalgraph - graph containing all edges
    %       Ggraphs - intermediate graphs
% NN = 1000; % number of edges, depends on how good the scores are.
% PCC_OVERLAP = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1
% 
% filtM = triu(PCC_OVERLAP);
% filtM(filtM==0) = nan;
% normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length
% 
% % normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
% [sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');

false_maps =[];

% barcodeGen = barcodeGen1;
G = digraph;
G = addnode(G,length(barcodeGen));
G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
Gundirected = graph;
Gundirected = addnode(Gundirected,length(barcodeGen));
Gundirected.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';

nodeOrientations = ones(1,length(barcodeGen));
 
Ggraphs =cell(1,NN);

edges =[];

for idxPair=1:NN
    % [maxScore,maxLoc] = max(normalizedScore);
    % idxPair = 1
    badBars = curStd>thrStd;
    if badBars(idxPair)~=1
        [curSink,curSource] = ind2sub(sizeOs,sortedIds(idxPair));
    
    %     import Core.plot_match_pcc;
    %     [f] = plot_match_pcc(barStruct, overlapStruct,curSink,curSource);
    
        % negative dist mean sink starts to the left of source. 
        dist = overlapStruct(curSink,curSource).pB;
%         orEdge = -2*overlapStruct(curSink,curSource).or+3;
        orEdge = overlapStruct(curSink,curSource).or; % if -1 1 already

        try
            path1 = shortestpath(Gundirected,curSource,curSink);
        catch
            path1 =[];
        end
%         try
%             path2 = shortestpath(G,curSink,curSource);
%         catch
%             path2 = [];
%         end
    
       pathExists = ~isempty(path1);%||~isempty(path2);
        if pathExists
            forwardMatch = (nodeOrientations(curSource)==orEdge)&&(nodeOrientations(curSink)==1);
            inverseMatch = (nodeOrientations(curSource)==(-1*orEdge))&&(nodeOrientations(curSink)==-1);
            matchingOrs = forwardMatch || inverseMatch;
%             if inverseMatch==1
%                 nodeOrientations(path1(2:end)) = -1*nodeOrientations(path1(2:end));
%             end
                
    %             matchingOrs
            % https://www.pnas.org/doi/pdf/10.1073/pnas.0604040103
            % elimination of false edges with inconsistent orientation
        else
            matchingOrs = 1;
            nearestNodes = nearest(Gundirected, curSink, Inf);
            if ~isempty(nearestNodes)
                nodeOrientations(nearestNodes) =  nodeOrientations(nearestNodes)*nodeOrientations(curSink)*orEdge*nodeOrientations(curSource);
            end
            nodeOrientations(curSink) = orEdge*nodeOrientations(curSource);  % and for all edges connected to curSink 
                        

        end
    
    
        if matchingOrs
            if dist > 0
                G = addedge(G, curSource, curSink, 1/dist);
            else
                G = addedge(G, curSource, curSink, -1/dist); % should be easy to show which one is source.
            end
            Gundirected = addedge(Gundirected, curSource, curSink);
            edges = [edges; curSource curSink];
        else
            false_maps = [false_maps; idxPair];
        end
    
    
            %
        Gtemp = G;
        
        listN = 1:length(barcodeGen);
        endNodes = Gtemp.Edges.EndNodes(:);
        curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';
        listN(curNodes) = [];
    
        Gtemp = rmnode(Gtemp,listN);
        weak_bins = conncomp(Gtemp,'Type','weak');
        numComp(idxPair) = max(weak_bins);
        numNodes(idxPair) = length(curNodes);
    
    
    % plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1)
    
    %     nexttile,plot(Gtemp,'Layout','force','ArrowSize',10,'MarkerSize',2)
    %   title(strcat([' numNodes = ' num2str(length(curNodes)) ', numComp = ' num2str(numComp(i))]))
      
      Ggraphs{idxPair} = Gtemp;
    end
end

nonempty = find(cellfun(@(x) ~isempty(x),Ggraphs));
if nargin < 8
    figure;
else
    h=tiledlayout(tsv,1,1);
    nexttile(h);
end
plot(Ggraphs{nonempty(end)},'Layout','force','ArrowSize',5,'MarkerSize',1)

finalgraph = Ggraphs{nonempty(end)};
end

