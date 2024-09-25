% Synth barcodes bargroup islands vs theory validation

sF = 0.9:0.025:1.1; % this should be the same as before?
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

nmPerPx = 110; % do we need?


%% recalculate positions for this single island / reference based
idx = 3;
bars = barcodeGen(barIslands{idx}); 

barx = barIslands{idx};

%todo: unique sf for each bar

sets.comparisonMethod = 'mass_pcc';

% compare to theory: either pcc or mp
[cSbi,rezMaxbi,~] = compare_to_t(bars,tS,sF,sets);

% comparison struct to overlap struct:
import Core.comparisonsres_to_struct; 
[oSreal] = comparisonsres_to_struct(bars,cSbi,tS{1}.length);


% synthS

% Plot: % don't have aligned kymo, so this plot is not applicable.. make a
% different plot?
% import CBT.Hca.UI.Helper.plot_any_bar;
% plot_any_bar(20,bars,rezMaxbi,tS,1);
% %

idxRun = 1
% bgen =  bG{idxRun};
[outConsensus, coverage, pval] = gen_reference_based_assembly(bars,cSbi,tS,'test11',inf);
% 
% usin GT
[outConsensus, coverage, pval] = gen_reference_based_assembly(bars,synthS(barIslands{idx}),tS,'test11',inf);
% comparison struct to overlap struct:
% import Core.comparisonsres_to_struct; 
% [oSreal2] = comparisonsres_to_struct(bars,synthS(barIslands{idx}),tS{1}.length);

%% CHECK #1: is there an overlap in the last 1000 pixels?
% put consensus barcodes into a structure
cGenEdge = [];
i = 1;

N = 1000;

cGenEdge{1}.rawBarcode = consensus{i}(1:N); % if longer than 2000
cGenEdge{1}.rawBitmask = ~isnan(consensus{i}(1:N));


cGenEdge{2}.rawBarcode = consensus{i}(end-N+1:end); % if longer than 2000
cGenEdge{2}.rawBitmask = ~isnan(consensus{i}(end-N+1:end));


w = 300;
% calculate MP scores 
[oSconsEdge{1}] = calc_overlap_mp(cGenEdge,sF, w,timestamp);
 

import Core.plot_match_islands;
[f] = plot_match_islands(cGenEdge, oSconsEdge{1}, [1:length(cGenEdge)-1],length(cGenEdge));



print('FIGS/FigSa.eps','-depsc','-r300');



%%

%

% NN = 2000;

bT{1}.rawBarcode = [bT{1}.rawBarcode];% bT{1}.rawBarcode(1:NN)]; %make a bit longer
% bT{1}.rawBitmask  = [bT{1}.rawBitmask bT{1}.rawBitmask(1:NN) ]

bT{1}.rawBitmask = ones(1,length(bT{1}.rawBarcode));
bT{1}.rawBitmask = logical(bT{1}.rawBitmask);

% put consensus barcodes into a structure
cGen = [];
edgePx = 1;
for i=1%:length(consensus)
    cGen{i}.rawBarcode = consensus{i}(1:end-edgePx+1);
    cGen{i}.rawBitmask = ~isnan(consensus{i}(1:end-edgePx+1));
end

numExpBar = length(consensus);

%  add theory as last barcode
cGen =[cGen bT];

w = 300;
% calculate MP scores 
[oScons{1}] = calc_overlap_mp(cGen,sF, w,timestamp);
 
% [oScons] = calc_scores({cGen},w,sF);

%%
% 
% import Core.plot_match_islands;
% [f] = plot_match_islands(cGen, oScons{1}, [1:length(cGen)-1],length(cGen));


import Core.plot_match_islands_oneplot;
[f,globconc] = plot_match_islands_oneplot(cGen,cGenAll, oScons{1}, [1:length(cGen)-1],length(cGen));

% [outConsensus, coverage, pval] = gen_reference_based_assembly(bars,cSbi,{theoryStruct},'test11',inf);

print('FIGS/FigSynth5.eps','-depsc','-r300');

%

%% Compare to theory results.
% import Core.comparisonsres_to_struct; 
% [rezStruct2] = comparisonsres_to_struct(barcodeGen, comparisonStructC{idxRun},thrS{idxRun}.length);
    
bIdx = 1;
posDifExp =  oScons{1}(bIdx,length(cGen)).pB -  oScons{1}(bIdx,length(cGen)).pA;

allPosLeft = cellfun(@(x) x.pos, cGenAll{bIdx}.comparisonStruct);
posStart = allPosLeft-min(allPosLeft);

posD = zeros(1,length(cGenAll{bIdx}.comparisonStruct));
for barNr = 1:length(cGenAll{bIdx}.comparisonStruct)
    ixCur = find(cGenAll{bIdx}.comparisonStruct{barNr}.barid==barx);
    posTHR =cSbi{ixCur}.pos(1) ;%+ comparisonStructC{idxRun}{cGenAll{bIdx}.comparisonStruct{barNr}.barid}.lengthMatch;
    
    posDifMatch  = posStart(barNr);
    if cSbi{ixCur}.or(1)==2
        posDifMatch  = length(cGenAll{1}.rawBarcode)-posDifMatch-length(bars{ixCur}.rawBarcode); % todo: add rescale for more accurate position
    end
    
    posD(barNr) = posTHR - (posDifExp + posDifMatch); % difference might be due to different length re-scaling factor w.r.t. theory
end

posD

f=figure
histogram(posD)
title('Pos. difference from mapping to ref. genome')
print('FIGS/Fig6.eps','-depsc','-r300');

super_quick_plot_rezmax(81,barcodeGen,rezMax,thrS{idxRun},1)

%% Plot specific pair
import Core.plot_match_simple;

expidx = 2; % idx of experiment
thrIdx = 1; % idx of theory
[f] = plot_match_simple(cGen, oScons{1},expidx,numExpBar+thrIdx);
[f] = plot_match_simple(cGen, oScons{1},numExpBar+thrIdx,expidx);

%
%% plot comparison for the best pair (this is both for theory vs exp and exp vs exp
idFold = 1;
idxPair = 1; % 1 is the pair with best score
 


import Core.plot_match_pcc;
[sortedVals, sortedIds,pscores] = sorted_scores(oScons{idFold});


[curSink,curSource] = ind2sub(size(oScons{idFold}),sortedIds(idxPair));

import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(cGen, oScons{idFold},curSink,curSource);


