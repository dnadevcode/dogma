% overlap graph construction 16/09/22
% modified 09/30/22

% PCC_OVERLAP = reshape([overlapStruct.score], size(overlapStruct,1),size(overlapStruct,2));
PCC_OVERLAP = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1

overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));
% lenOverlap = 

%%
%% filter scores based on best overlap
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');

% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 30
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,xId, yId);

%%
NN = 1000; % number of edges, depends on how good the scores are.



barcodeGen = barcodeGen1;
G = digraph;
G = addnode(G,length(barcodeGen));
G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';

nodeOrientations = ones(1,length(barcodeGen));
 
Ggraphs =cell(1,NN);

edges =[];

for idxPair=1:NN
    % [maxScore,maxLoc] = max(normalizedScore);
    % idxPair = 1
    badBars = curStd>thrStd;
    if badBars(idxPair)~=1
        [curSink,curSource] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    
    %     import Core.plot_match_pcc;
    %     [f] = plot_match_pcc(barStruct, overlapStruct,curSink,curSource);
    
        % negative dist mean sink starts to the left of source. 
        dist = overlapStruct(curSink,curSource).pB;
        orEdge = -2*overlapStruct(curSink,curSource).or+3;
    
        try
            path1 = shortestpath(G,curSource,curSink);
        catch
            path1 =[];
        end
        try
            path2 = shortestpath(G,curSink,curSource);
        catch
            path2 = [];
        end
    
       pathExists = ~isempty(path1)||~isempty(path2);
        if pathExists
            matchingOrs = (nodeOrientations(curSource)==orEdge)&&(nodeOrientations(curSink)==1)||...
                (nodeOrientations(curSource)==(-1*orEdge))&&(nodeOrientations(curSink)==-1);
    %             matchingOrs
            % https://www.pnas.org/doi/pdf/10.1073/pnas.0604040103
            % elimination of false edges with inconsistent orientation
        else
            matchingOrs = 1;
            nodeOrientations(curSink) = orEdge;  
    
        end
    
    
        if matchingOrs
            if dist > 0
                G = addedge(G, curSource, curSink, 1/dist);
            else
                G = addedge(G, curSource, curSink, -1/dist); % should be easy to show which one is source.
            end
            edges = [edges; curSource curSink];
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
figure
plot(Ggraphs{nonempty(end)},'Layout','force','ArrowSize',5,'MarkerSize',1)
%%

figure
plot(G,'Layout','force','ArrowSize',5,'MarkerSize',1)

% [bins,binsizes] = conncomp(Ggraphs{end});

[bins,binsizes] = conncomp(G,'Type','weak');
idx = binsizes(bins) == max(binsizes);
SG = subgraph(G,idx);

nodes = cellfun(@(x) str2num(x),SG.Nodes.Variables);

nonemptyrows = find(arrayfun(@(x) ~isempty(intersect(edges(x,:),nodes)),1:size(edges,1)));


figure
plot(SG,'Layout','force','ArrowSize',5,'MarkerSize',1)


listOfBarPairsInBargroup = edges(nonemptyrows,:)


foer ii=1:length(nonemptyrows)
fig =figure;

    ix2 = edges(nonemptyrows(ii),1);
    ix1 = edges(nonemptyrows(ii),2);
    plot([-overlapStruct(ix1,ix2).pB+2 -overlapStruct(ix1,ix2).pB+2+overlapStruct(ix1,ix2).lenB-1],[2*i 2*i],'black|-')
    hold on
    plot([overlapStruct(ix1,ix2).pA overlapStruct(ix1,ix2).pA+overlapStruct(ix1,ix2).lenA-1],[2*i+1 2*i+1],'black|-')
ylim([0 4])

%% now we plot overlap plot. For two barcodes, it's just their relative position


fig =figure;
ix1=141;
ix2 = 73;
ix3 = 131;
plot([overlapStruct(ix1,ix2).pB overlapStruct(ix1,ix2).pB+overlapStruct(ix1,ix2).lenB-1],[1 1],'black|-')
hold on
plot([overlapStruct(ix1,ix2).pA overlapStruct(ix1,ix2).pA+overlapStruct(ix1,ix2).lenA-1],[2 2],'black|-')
% plot([overlapStruct(ix3,ix2).pA overlapStruct(ix3,ix2).pA+overlapStruct(ix3,ix2).lenA-1],[3 3],'black|-')

ylim([0 4])

import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,ix1, ix2);


% 
%     % check if edge already exists (since for now we are comparing both
%     % A->B and B->A
% %     [idxOut,m] = findedge(G,curSource,curSink);
% %     [idxOut2,m] = findedge(G,curSink,curSource);
%    
% 
% 
%     % add edge
%     if idxOut==0 && idxOut2==0     
%         % orientation consistency check. If path between source and sink
%         % exists, the two nodes belong to the same component of the graph.
%         % Orientation of maps represented by these nodes must be consistent
%         % with orientation of maps with suggested overlap
%         try
%             path1 = shortestpath(G,curSource,curSink);
%         catch
%             path1 =[];
%         end
%         try
%         path2 = shortestpath(G,curSink,curSource);
%         catch
%             path2 = [];
%         end
% 
%         pathExists = ~isempty(path1)||~isempty(path2);
%         if pathExists
%             matchingOrs = (nodeOrientations(curSource)==barOrs{i}(1))&&(nodeOrientations(curSink)==1)||...
%                 (nodeOrientations(curSource)==(-1*barOrs{i}(1)))&&(nodeOrientations(curSink)==-1);
% %             matchingOrs
%             % https://www.pnas.org/doi/pdf/10.1073/pnas.0604040103
%             % elimination of false edges with inconsistent orientation
%         else
%             matchingOrs = 1;
%             nodeOrientations(curSource) = barOrs{i}(1);  
% 
%         end
%         
%         if matchingOrs
%             if dist > 0
%                 G = addedge(G, curSource, curSink, 1/dist);
%             else
%                 G = addedge(G, curSink, curSource, -1/dist);
%             end
%         end
%     end
% 
% 
% %
% 
% G = digraph;
% G = addnode(G,length(barcodeGen));
% nodeOrientations = ones(1,length(barcodeGen));
%   
% for i=1:length(scoresGood)
%     curSource = barsSource(scoresIdx(i));
%     curSink =  barsSink(scoresIdx(i));
%     
%     import Core.plot_bargroup;
%     [barMat{i},bars{i},barOrs{i}] =  plot_bargroup(curSource,stridx,mpI1(curSink),LOCS1(curSink),pksUnique1(curSink),...
%         pksUniquePos1(curSink),curSink,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
% 
%     center1 = mean(barMat{i}{1}{1});
%     center2 = mean(barMat{i}{2}{1});
% %     barOrs{i}(j)
%     dist = center1 - center2;
%     
%     % check if edge already exists (since for now we are comparing both
%     % A->B and B->A
%     [idxOut,m] = findedge(G,curSource,curSink);
%     [idxOut2,m] = findedge(G,curSink,curSource);
%         
%      
%     % add edge
%     if idxOut==0 && idxOut2==0
%         
%         % orientation consistency check. If path between source and sink
%         % exists, the two nodes belong to the same component of the graph.
%         % Orientation of maps represented by these nodes must be consistent
%         % with orientation of maps with suggested overlap
%         try
%         path1 = shortestpath(G,curSource,curSink);
%         catch
%             path1 =[];
%         end
%         try
%         path2 = shortestpath(G,curSink,curSource);
%         catch
%             path2 = [];
%         end
% 
%         pathExists = ~isempty(path1)||~isempty(path2);
%         if pathExists
%             matchingOrs = (nodeOrientations(curSource)==barOrs{i}(1))&&(nodeOrientations(curSink)==1)||...
%                 (nodeOrientations(curSource)==(-1*barOrs{i}(1)))&&(nodeOrientations(curSink)==-1);
% %             matchingOrs
%             % https://www.pnas.org/doi/pdf/10.1073/pnas.0604040103
%             % elimination of false edges with inconsistent orientation
%         else
%             matchingOrs = 1;
%             nodeOrientations(curSource) = barOrs{i}(1);  
% 
%         end
%         
%         if matchingOrs
%             if dist > 0
%                 G = addedge(G, curSource, curSink, 1/dist);
%             else
%                 G = addedge(G, curSink, curSource, -1/dist);
%             end
%         end
%     end
% end


f=figure,plot(G,'Layout','force','ArrowSize',5,'MarkerSize',2)
saveas(f,'figs/fig5.png')


