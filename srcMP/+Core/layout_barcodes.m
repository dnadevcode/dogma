function [outputArg1,outputArg2] = layout_barcodes(inputArg1,inputArg2)

% this is the layout step of bargrouping prototype method


% we have a p-val thresh
pthresh = 0.0001;
%%
pvals = cell(1,length(barcodeGen));

% G = graph;
% G = addnode(G,length(barcodeGen));
% % 
% for i=1:length(barcodeGen)
% end

% calculate p-vals

for i=1:length(barcodeGen)
    i
    bar1 = i;
    
    % pcc values for this barcode
    allPCCvals = PKS1{bar1}(pksUniquePos1{bar1});
    allPCCvals=allPCCvals(allPCCvals>0.75); % reduce number of things we map

    % PVAL for thes pcc values
    len2 = lengths(pksUnique1{bar1});
    len1 = lengths(bar1);
    
    % pval params. for specific PSF and nm/px ratios//
    % pval_test_loop/pval_test_real 
    pvalPar1 = 0.13*MIN_OVERLAP_PIXELS;
    pvalPar2 = 2*(len1-MIN_OVERLAP_PIXELS)*(len2-MIN_OVERLAP_PIXELS)/100;
%     pval =
    import Zeromodel.beta_ev_pdf;

%     [p] = beta_ev_pdf(allPCCvals', pvalPar1, 1, pvalPar2);%y_est(k)
    p = ones(1,length(allPCCvals));
    for j=1:length(allPCCvals)
        [p(j)] =  1-Zeromodel.beta_ev_cdf(allPCCvals(j), pvalPar1, 1, pvalPar2(j), 1);
%         if p(j) > 0.5 % early stopping criteria, sometimes not right? 
%             break;
%         end
    end
    
    pthresh = 0.00001;
    numGoodMatch = p < pthresh;
    peaksToTry = pksUnique1{bar1}(numGoodMatch);
    pvals{i} = p;
    pleaksI{i} = peaksToTry;

%    graphMatrix(bar1,peaksToTry) = 1; % just the graph matrix
%     graphMatrix(peaksToTry,bar1) = 1; % just the graph matrix

%     allb(peaksToTry)
    % now create bargroup for bar1.
    import Core.plot_bargroup;
    [barMat{i},bars{i},barOrs{i}] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
     %                                       ix,stridx,mpI                   ,LOCS                   ,pksUnique1                  ,pksUniquePos                   ,k,               baridx, sF, barcodeGenGood, h
%     G = graph;
%      G = addnode(G,i);

end

%% create overlap graph

% First sort the overlaps based on the score
pthresh = 0.0001;

scoresGood = [];
barsSource = [];
barsSink = [];
barsPCC = [];
for i=1:length(barcodeGen)
    bar1 = i;
    allPCCvals = PKS1{bar1}(pksUniquePos1{bar1});
    allPCCvals=allPCCvals(allPCCvals>0.75); % reduce number of things we map

    numGoodMatch =  pvals{i} < pthresh;
    peaksToTry = pksUnique1{bar1}(numGoodMatch);

    barsPCC = [barsPCC; allPCCvals(numGoodMatch)];
    scoresGood = [scoresGood pvals{i}(pvals{i} < pthresh)];
    barsSource = [ barsSource i*ones(1,length(peaksToTry))];
    barsSink = [ barsSink peaksToTry];
end

[a,scoresIdx] = sort(scoresGood);

% [a,scoresIdx] = sort(barsPCC,'desc');
% barsPCC(scoresIdx)

%% now main loop through the scores..

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
        path1 = shortestpath(G,curSource,curSink);
        path2 = shortestpath(G,curSink,curSource);

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
                G = addedge(G, curSource, curSink, dist);
            else
                G = addedge(G, curSink, curSource, -dist);
            end
        end
    end
end

figure,plot(G,'Layout','force','ArrowSize',5,'MarkerSize',2)

% figure,plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'ArrowSize',5,'MarkerSize',2)
% figure,plot(G)
% 
% [bin,binsize] = conncomp(G,'Type','weak');
% idx = binsize(bin) == max(binsize);
% SG = subgraph(G,idx);
% figure;plot(SG)
% Elimination of False Edges with Consistent Orientation..

%%
G = digraph;
  G = addnode(G,2*length(barcodeGen));
% 
% for i=1:length(barcodeGen)
% end

for i=1:length(barcodeGen)
        bar1 = i;
    numGoodMatch =  pvals{i} < pthresh;
    peaksToTry = pksUnique1{bar1}(numGoodMatch);
    import Core.plot_bargroup;
    [barMat{i},bars{i},barOrs{i}] =  plot_bargroup(i,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);

    for j=1:length(peaksToTry)
%         barMat{i}{j}
      
        import Core.plot_bargroup;
        plot_bargroup(i,stridx,mpI1(barMat{i}{j+1}{3}),LOCS1(barMat{i}{j+1}{3}),pksUnique1(barMat{i}{j+1}{3}), pksUniquePos1(barMat{i}{j+1}{3}),barMat{i}{j+1}{3},baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,1);
        
        center1 = mean(barMat{i}{1}{1});
        center2 = mean(barMat{i}{j+1}{1});
        barOrs{i}(j)
        dist = center1 - center2;
        
        % we want to add edge between i,or_i -> j if dist>0, and
        % j -> i,or_i if dist < 0
        
        % note here that j has a fixed orientation, and i can be forward or
        % reverse.
        
        %
        
        
        
        
        if dist > 0

            [idxOut,m] = findedge(G,i,peaksToTry(j));
            if idxOut==0
                G = addedge(G,i,peaksToTry(j),barOrs{i}(j));
            end
        end
    end
end

%%
f=figure;
plot(G)
% H = transreduction(G)
% f=figure
% plot(H)

%%

G = digraph;
  G = addnode(G,length(barcodeGen));
% 
% for i=1:length(barcodeGen)
% end

for i=1:length(barcodeGen)
    
    for j=1:length(pleaksI{i})
%         G = addnode(G,peaksToTry(j));
%         G.
    [idxOut,m] = findedge(G,i,pleaksI{i}(j));
    [idxOut2,m] = findedge(G,pleaksI{i}(j),i);

    if idxOut==0 && idxOut2==0
        if barOrs{i}(j) == 1;
            G = addedge(G,i,pleaksI{i}(j));
        else
            G = addedge(G, pleaksI{i}(j), i);
        end
    end
end
end

f=figure
plot(G)

H = transreduction(G);

figure,plot(H)

%%
% peaksToTry
barIdx = peaksToTry;

% if theory generated..
[outConsensus] = gen_reference_based_assembly(barcodeGen(barIdx),comparisonStruct(barIdx),theoryStruct);

 %%
allbars = b;
[c,idPos] = sort(b);
allb = ones(1,length(allbars));
%
bvec= [];

graphMatrix = zeros(length(barcodeGen),length(barcodeGen));

for ii=1:10
    if sum(allb) ==0
        break;
    end
    % now loop through first few
    btemp = allbars(find(allb));
    idx = 1;
    bar1 = btemp(idx);

    % take PCC vals for unique matching barcodes
    allPCCvals = PKS1{bar1}(pksUniquePos1{bar1}); %% CHECK/POSSIBLE BUG
    
    % convert to PVAL
    len2 = lengths(pksUnique1{bar1});
    len1 = lengths(bar1)
    pvalPar1 = 0.13*MIN_OVERLAP_PIXELS;
    pvalPar2 = 2*(len1-MIN_OVERLAP_PIXELS)*(len2-MIN_OVERLAP_PIXELS)/100;
    import Zeromodel.beta_ev_pdf;
%     pval =

%     [p] = beta_ev_pdf(allPCCvals', pvalPar1, 1, pvalPar2);%y_est(k)
    for j=1:length(allPCCvals)
        [p(j)] =  1-Zeromodel.beta_ev_cdf(allPCCvals(j), pvalPar1, 1, pvalPar2(j), 1);
    end

    
    % find number above thresh
%     numGoodMatch = allPCCvals > 0.8; % here for PCC
    numGoodMatch = p < 0.01;
    % threshold peaks based on intensity 
    % % now for a single barcode check the cluster barcodes
    % bar1= goodIdx(8);  % 709
    % bar1= 548;  % 709
    peaksToTry = pksUnique1{bar1}(numGoodMatch);
    
    graphMatrix(bar1,peaksToTry) = 1; % just the graph matrix
%     graphMatrix(peaksToTry,bar1) = 1; % just the graph matrix

%     allb(peaksToTry)
    % now create bargroup for bar1.
%     import Core.plot_bargroup;
%     [barMat{ii},bars{ii}] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
     %                                       ix,stridx,mpI                   ,LOCS                   ,pksUnique1                  ,pksUniquePos                   ,k,               baridx, sF, barcodeGenGood, h

    allb(idPos([bar1 peaksToTry])) = 0;
%     bvec{ii} = [bar1 peaksToTry];
end
    figure,imagesc(graphMatrix)
    
    
end

