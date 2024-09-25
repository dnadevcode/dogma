function [outputArg1,outputArg2] = layout_consensus(G)
% from overlap graph we create a layout

figure,plot(G,'Layout','force','ArrowSize',5,'MarkerSize',2,'NodeColor','red')


% firstSource = barsSource(1);
% firstSink = barsSink(1);

% allpaths(G,firstSink)

% If graph is a circle, we must find a cycle through the subgraph
% maximizing the genomic distance spanned by this path

% first create subraph
[bins,binsizes] = conncomp(G);
idx = binsizes(bins) == max(binsizes);
SG = subgraph(G,idx);

% nodesInMainComponent = find(idx); % alternatively could be a few components

component_sources = find(idx);
path_length=[];path_idx=[];
d= [];
% v = cell(1,length(component_sources));
for i=1:length(component_sources)
    i
%     v{i} = dfsearch(G,firstSource,'allevents');
    % find shortest path three, alternative to depth-first search
    [TR,D] = shortestpathtree(G,component_sources(i));
    % restore the weights, since for shortest path we used 1/distance
    TR.Edges.Weight=1./TR.Edges.Weight;
    d{i}=distances(TR);
    d{i}(isinf( d{i}))=0;

    [path_length(i),path_idx(i)]=max( d{i}(:));
end

[a,idx_a] = max(path_length);

[b,idx_b] = max(d{idx_a}(:));

[idx,idy] = ind2sub(size(d{idx_a}),idx_b);
path = shortestpath(G,idx,idy);

f=figure,plot(TR,'Layout','force','ArrowSize',5,'MarkerSize',2,'NodeColor','red')
save('figs/fig6.png')
% now we can use the 
%%
draftMap = [];
for i=1:length(path)-1
curSource = path(i);
curSink = path(i+1);
   import Core.plot_bargroup;
    [barMatdraft{i},barsdraft{i},barOrsdraft{i},reFacdraft{i}] =  plot_bargroup(curSource,stridx,mpI1(curSink),LOCS1(curSink),pksUnique1(curSink),...
        pksUniquePos1(curSink),curSink,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
    draftMap =[ draftMap barMatdraft{i}{2}{2}(barMatdraft{i}{2}{1}<0)];
end
draftMap = [draftMap barMatdraft{i}{1}{2}(barMatdraft{i}{2}{1}>0)];


% [bar1idx,bar2idx] = [29],[49];
vec = [29 65];
ix=5;
vec = path(ix:ix+1);

% vec= [248 29];
bar1idx=vec(1);
bar2idx=vec(2);
sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);
import Core.plot_best_match;
[f,pos1, bar1, pos2, bar2,pcc1,pcc2] = plot_best_match(sortedidx2,stridx,mpI1{bar1idx},LOCS1{bar1idx},pksUniquePos1{bar1idx},bar1idx,baridx2{bar1idx},sF,barStruct,MIN_OVERLAP_PIXELS);


barcodeGenDraft{1}.rawBarcode =  draftMap;
barcodeGenDraft{1}.rawBitmask =  ones(1,length(draftMap));


% todo: compare in the same way as MP. 
% can create one long file from theoryStruct.
rezMaxDraft=[];bestBarStretchDraft=[];bestLengthDraft=[];
for i=1:length(theoryStruct)
 tic
 import CBT.Hca.Core.Comparison.on_compare;
[rezMaxDraft{i},bestBarStretcHdRAFT{i},bestLengthdRAFT{i}] = on_compare(barcodeGenDraft,theoryStruct{i},'mass_pcc',sF,[],50);
toc
end

% MAP refinement procedure


% 
% figure,p=plot(SG)
% highlight(p,TR,'EdgeColor','r')

% depth-first search

end

