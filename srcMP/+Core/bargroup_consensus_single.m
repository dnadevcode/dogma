

% should we have two graphs? 
% one graph for bargroups
% second graph for all barcodes

% so bargroup graph links the relevant bargroups, where the second graph is
% the one used to check of two bargroups are overlapping. We merge (either that or discard one of them)
% two bargroups together if they overlap (cover the same area) and are not
% differing significantly. 


% INITIALIZATION

filtM = triu(PCC_OVERLAP);

filtM;
filtM(filtM==0) = nan;

pthresh = 0.00001;


G = digraph;
G = addnode(G,length(barcodeGen));
G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
nodeOrientations = ones(1,length(barcodeGen));

N = 100;

B = digraph;
B = addnode(B,N);

% f = figure
% tiledlayout(ceil(sqrt(N)),ceil(N/ceil(sqrt(N))),'TileSpacing','compact','Padding','compact')

Ggraphs = cell(1,N);

nodes = cell(1,N);
barMvec =  cell(1,N);
for i=1:N
    i
    bargroupMerge = 0;
    % [maxScore,maxLoc] = max(normalizedScore);
    idxPair = 1;

    normalizedScore = filtM(:).*sqrt(lenOverlap(:));
    [sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');


    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

    [PCC_OVERLAP(xId,yId) xId yId]
    bar1idx=xId;bar2idx=yId;
    sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);
        import Core.plot_best_match;
        [f,pos1, bar1, pos2, bar2,pcc1,pcc2] = plot_best_match(sortedidx2,stridx,mpI1{bar1idx},LOCS1{bar1idx},pksUniquePos1{bar1idx},bar1idx,baridx2{bar1idx},sF,barStruct,MIN_OVERLAP_PIXELS);
        pcc1
    % pcc2

    % xId = 104
%     import Core.single_bargroup;
%      [barMat,bars,orBars,reFac,pval] = single_bargroup(xId,pthresh,PCC_OVERLAP,PCC_MP,len1,len2,lenOverlap,...
%             stridx,mpI1,LOCS1,pksUnique1,pksUniquePos1,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,1);
    

    nodes{i} = cellfun(@(x) str2num(x),bars);
    % bargroup map
    import Core.bargroup_map;
     [barMat,shiftVec,barM,barMzscored] = bargroup_map(bars, reFac, orBars,barcodeGen,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers);
    barMvec{i} = barM;
    barMveczcored{i} = barMzscored;


    % check if bargroup overlaps with previous bargroup
     for j=1:i-1
         if ~isempty(intersect(nodes{i},nodes{j}))
            % then add edge, unless containment is inconsistent 
            
            %check overlap between bargroups
            B = addedge(B, i, j, 1);
%             
%             bargroup1 = barMvec{j}
%             bargroup2 = barMvec{i}
%    ;
    bargroup1 = barMveczcored{47};
    bargroup2 = barMveczcored{75};

    bargroupGen{1}.rawBarcode = nanmean(bargroup1,1);
    bargroupGen{1}.rawBitmask = sum(~isnan(bargroup1(:,1:end-1)))>=2;
    bargroupGen{1}.rawBitmask(end+1) = 0;

    bargroupGen{2}.rawBarcode = nanmean(bargroup2,1);
    bargroupGen{2}.rawBitmask = sum(~isnan(bargroup2(:,1:end-1)))>=2;
    bargroupGen{2}.rawBitmask(end+1) = 0;
       
    import Core.bargroup_map;
    [barMat,shiftVec,barM,barMzscored] = bargroup_map({'1','2'}, 1, 1,bargroupGen,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers);
  
     figure,plot(barM')
     figure,plot(barMzscored')

%              figure,plot()

%                 import Core.bargroup_map;
%      [barMat,shiftVec,barM,barMzscored] = bargroup_map(bars, reFac, orBars,barcodeGen,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers);
 
         end
     end
     
       Ggraphs{i} = B;

     
    path = cellfun(@(x) str2num(x),bars);
    allpairs = nchoosek(path,2);

    for tt=1:size(allpairs,1)
        filtM(allpairs(tt,1),allpairs(tt,2)) = nan;
        filtM(allpairs(tt,2),allpairs(tt,1)) = nan;
    end
     
end
%%
figure
plot(Ggraphs{100},'Layout','force','ArrowSize',5,'MarkerSize',1)

         %
Gtemp = B;

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
  
% coverage for this bargroup.
% figure,plot(sum(~isnan(barM(:,1:end-1))))

% remove all pairs with


% end
%     % center-1
%     center1 = mean(barMat{1}{1});
%     curSource = xId;
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
%         try %only need one way with unweifhted
%             path1 = shortestpath(G,curSource,curSink,'method','unweighted');
%         catch
%             path1 =[];
%         end
%         try
%             path2 = shortestpath(G,curSink,curSource,'method','unweighted');
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
%             bargroupMerge = 1;
%         else
            matchingOrs = 1;
            nodeOrientations(curSource) = orBars(j-1);  
            if bargroupMerge == 0
                B = addnode(B,1);
            end
            bargroupMerge = 1;

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
  
% coverage for this bargroup.
% figure,plot(sum(~isnan(barM(:,1:end-1))))

% remove all pairs with
    path = cellfun(@(x) str2num(x),bars);
    allpairs = nchoosek(path,2);

    for tt=1:size(allpairs,1)
        filtM(allpairs(tt,1),allpairs(tt,2)) = nan;
        filtM(allpairs(tt,2),allpairs(tt,1)) = nan;
    end
% filtM([allpairs(:,2) allpairs(:,1)]) = nan;
% end