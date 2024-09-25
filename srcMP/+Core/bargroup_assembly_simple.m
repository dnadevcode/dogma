function [Ggraphs,numComp,numNodes] = bargroup_assembly_simple(PCC_OVERLAP,PCC_MP, lenOverlap, len1,len2,pksUnique1,sortedidx2,stridx,mpI1,LOCS1,...
    pksUniquePos1,bar1idx,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,barcodeGen,pthresh)
    % bargroup based asembly.
    
    import Core.single_bargroup;

    %%
    PCC_OVERLAP = reshape([overlapStruct.fullscore], size(overlapStruct,1),size(overlapStruct,2))
    overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2))
    % lenOverlap = 
    
    %%
    %% filter scores based on best overlap
    filtM = triu(PCC_OVERLAP);
    filtM(filtM==0) = nan;
    normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
    [sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');
% 
% % [maxScore,maxLoc] = max(normalizedScore);
% idxPair = 1;
% [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
% 
% 
% PCC_OVERLAP(xId,yId)

% import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, overlapStruct,xId,yId);
% 
%  
% pthresh = 0.00001;
% import Core.simple_bargroup;
% [orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,xId,pthresh,1);

    
G = digraph;
G = addnode(G,length(barcodeGen));
G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
graphMatrix = zeros(length(barcodeGen),length(barcodeGen));

nodeOrientations = ones(1,length(barcodeGen));
  
N = 64;

% f = figure
% tiledlayout(ceil(sqrt(N)),ceil(N/ceil(sqrt(N))),'TileSpacing','compact','Padding','compact')

[allxId,allyId] = ind2sub(size(PCC_OVERLAP),sortedIds);

Ggraphs = cell(1,N);
for i=1:N
%     idxPair = i;
    xId = allxId(1);
%     yId = allyId(1);
%     [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

    curSource = xId;

    
    pthresh = 0.00001;
    import Core.simple_bargroup;
    [orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,curSource,pthresh,1);
    [res] = get_display_results(barcodeGen(bars),[], comparisonStruct(bars), theoryStruct, sets);
    import CBT.Hca.UI.Helper.plot_best_pos;
    plot_best_pos([], compStr(bars), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

    graphMatrix(curSource,bars(2:end)) = 1; % just the graph matrix\

    for k=1:length(bars)
        allyId(allxId==bars(k))=[];
        allxId(allxId==bars(k))=[];
    end


    for j=1:length(bars)-1
        G = addedge(G, curSource, bars(j+1), 1);
    end
    
            %
    Gtemp = G;
    
    listN = 1:length(barcodeGen);
    endNodes = Gtemp.Edges.EndNodes(:);
    curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';
    listN(curNodes) = [];

    Gtemp = rmnode(Gtemp,listN);
    weak_bins = conncomp(Gtemp,'Type','weak');
    figure, plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1)

%     numComp(i) = max(weak_bins);
%     numNodes(i) = length(curNodes);

% 
%     % center-1
%     center1 = mean(barMat{1}{1});
% 
%     for j=2:length(barMat)
%         curSink =  barMat{j}{3};
%     
%         % center-2
%         center2 = mean(barMat{j}{1});
%         
%         % distance between centers
%         dist = center1 - center2;
%     
%         % check if edge already exists (since for now we are comparing both
%         % A->B and B->A
%         [idxOut,m] = findedge(G,curSource,curSink);
%         [idxOut2,m] = findedge(G,curSink,curSource);
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
%             path1 = shortestpath(G,curSource,curSink);
%         catch
%             path1 =[];
%         end
%         try
%             path2 = shortestpath(G,curSink,curSource);
%         catch
%             path2 = [];
%         end
% 
%         pathExists = ~isempty(path1)||~isempty(path2);
%         if pathExists
%             matchingOrs = (nodeOrientations(curSource)==orBars(j-1))&&(nodeOrientations(curSink)==1)||...
%                 (nodeOrientations(curSource)==(-1*orBars(j-1)))&&(nodeOrientations(curSink)==-1);
% %             matchingOrs
%             % https://www.pnas.org/doi/pdf/10.1073/pnas.0604040103
%             % elimination of false edges with inconsistent orientation
%         else
%             matchingOrs = 1;
%             nodeOrientations(curSource) = orBars(j-1);  
%         end
%         
%         if matchingOrs % add positive dist
%             if dist > 0
%                 G = addedge(G, curSource, curSink, 1/dist);
%             else
%                 G = addedge(G, curSink, curSource, -1/dist);
%             end
%         end
%     end
%     end

    %
Gtemp = G;

listN = 1:length(barcodeGen);
endNodes = Gtemp.Edges.EndNodes(:);
curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';
listN(curNodes) = [];

    Gtemp = rmnode(Gtemp,listN);
    weak_bins = conncomp(Gtemp,'Type','weak');
    numComp(i) = max(weak_bins);
    numNodes(i) = length(curNodes);

% plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1)

%     nexttile,plot(Gtemp,'Layout','force','ArrowSize',10,'MarkerSize',2)
%   title(strcat([' numNodes = ' num2str(length(curNodes)) ', numComp = ' num2str(numComp(i))]))
  
  Ggraphs{i} = Gtemp;
  
end
%%
% saveas(f,'figs/figOverlapGraph.png')

%%
Gtemp = G;

listN = 1:length(barcodeGen);
endNodes = Gtemp.Edges.EndNodes(:);
curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';

listN(curNodes) = [];

    Gtemp = rmnode(Gtemp,listN);

figure;
plot(Gtemp,'Layout','force','ArrowSize',10,'MarkerSize',1)
% 
% figure
% P = plot(Gtemp,'LineWidth',Gtemp.Edges.Weight);
