function [pval,G] = layout_graph(MIN_OVERLAP_PIXELS,PCC_OVERLAP,len1,len2)

% [PCC_OVERLAP,len1,len2
%         min_overlap = lengthsbar(i); %min(lengthsThry(i),lengthsbar(i));

    a = 0.13*MIN_OVERLAP_PIXELS;

    pval = ones(size(len1));

% first gen p-values, based on the full overlap (alternative, more prone to misalignment: based on
% MIN_OVERLAP)

    for i=1:size(len1,1)
i
        for j=1:size(len1,2)
            n = 2*(len2(i,j)+len1(i,j)-2*MIN_OVERLAP_PIXELS);
            if PCC_OVERLAP(i,j) > 0
                pval(i,j) = 1-Zeromodel.beta_ev_cdf(PCC_OVERLAP(i,j), a, 1, n, 1);
            end
        end
    end

%%
pthresh = 0.1;

scoresGood = [];
barsSource = [];
barsSink = [];
barsPCC = [];

for i=1:length(len1)
% i
    % find good matches
    pvalCur = pval(i,:);
    numGoodMatch = find(pval(i,:)<pthresh);

    peaksToTry = numGoodMatch;



    barsPCC = [barsPCC; PCC_OVERLAP(i,numGoodMatch)'];
    scoresGood = [scoresGood pvalCur(numGoodMatch)];
    barsSource = [ barsSource i*ones(1,length(peaksToTry))];
    barsSink = [ barsSink peaksToTry];
end

    [a,scoresIdx] = sort(scoresGood);


% now main loop through the scores..

G = digraph;
G = addnode(G,length(barcodeGen));
nodeOrientations = ones(1,length(barcodeGen));
  
for i=1:length(scoresGood)
    curSource = barsSource(scoresIdx(i));
    curSink =  barsSink(scoresIdx(i));
    
    import Core.plot_bargroup;
    [barMat{i},bars{i},barOrs{i}] =  plot_bargroup(curSource,stridx,mpI1(curSink),LOCS1(curSink),pksUnique1(curSink),...
        pksUniquePos1(curSink),curSink,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);

    center1 = mean(barMat{i}{1}{1});
    center2 = mean(barMat{i}{2}{1});
%     barOrs{i}(j)
    dist = center1 - center2;
    
    % check if edge already exists (since for now we are comparing both
    % A->B and B->A
    [idxOut,m] = findedge(G,curSource,curSink);
    [idxOut2,m] = findedge(G,curSink,curSource);
        
     
    % add edge
    if idxOut==0 && idxOut2==0
        
        % orientation consistency check. If path between source and sink
        % exists, the two nodes belong to the same component of the graph.
        % Orientation of maps represented by these nodes must be consistent
        % with orientation of maps with suggested overlap
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
            matchingOrs = (nodeOrientations(curSource)==barOrs{i}(1))&&(nodeOrientations(curSink)==1)||...
                (nodeOrientations(curSource)==(-1*barOrs{i}(1)))&&(nodeOrientations(curSink)==-1);
%             matchingOrs
            % https://www.pnas.org/doi/pdf/10.1073/pnas.0604040103
            % elimination of false edges with inconsistent orientation
        else
            matchingOrs = 1;
            nodeOrientations(curSource) = barOrs{i}(1);  

        end
        
        if matchingOrs
            if dist > 0
                G = addedge(G, curSource, curSink, 1/dist);
            else
                G = addedge(G, curSink, curSource, -1/dist);
            end
        end
    end
end

f=figure,plot(G,'Layout','force','ArrowSize',5,'MarkerSize',2)
saveas(f,'figs/fig5.png')




end

