
% barislands vs theory when each barcode is considered individually.


w = 300;% minimum overlap length
sF = 0.9:0.025:1.1;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');


% quickly calculate theory for given single fasta file
[theoryStructRev,theoryStruct,bT] = prep_thry_for_local_comp(fastaFile, nmbp, nmPerPx, 1);

bT{1}.rawBitmask = logical(bT{1}.rawBitmask);

bars = outConsensus{2};

% put consensus barcodes into a structure
cGen = [];
for i=1:size(bars,1)
    cGen{i}.rawBarcode = bars(i,:);
    cGen{i}.rawBitmask = ~isnan( bars(i,:));
end

% numExpBar = length(consensus);

%  add theory as last barcode
% cGen =[cGen bT];

w = 500;
% calculate MP scores 
[oScons{1}] = calc_overlap_mp_multidim(cGen,bT,sF, w,timestamp);
 
% [oScons] = calc_scores({cGen},w,sF);

%%

import Core.plot_match_islands;
[f] = plot_match_islands(cGen, oScons{1}, [1:length(cGen)-1],length(cGen));


import Core.plot_match_islands_oneplot;
[f,globconc] = plot_match_islands_oneplot(cGen,cGenAll, oScons{1}, [1:length(cGen)-1],length(cGen));

[outConsensus, coverage, pval] = gen_reference_based_assembly(bgen,comparisonStructC{idxRun},thrS(idxRun),'test11',inf);

print('FIGS/Fig5.eps','-depsc','-r300');

%% Compare to theory results.
% import Core.comparisonsres_to_struct; 
% [rezStruct2] = comparisonsres_to_struct(barcodeGen, comparisonStructC{idxRun},thrS{idxRun}.length);
    
bIdx = 7;
posDifExp =  oScons{1}(bIdx,length(cGen)).pB -  oScons{1}(bIdx,length(cGen)).pA;

allPosLeft = cellfun(@(x) x.pos, cGenAll{bIdx}.comparisonStruct);
posStart = allPosLeft-min(allPosLeft);

posD = zeros(1,length(cGenAll{bIdx}.comparisonStruct));
for barNr = 1:length(cGenAll{bIdx}.comparisonStruct)
    ixCur = cGenAll{bIdx}.comparisonStruct{barNr}.barid;
    posTHR = comparisonStructC{idxRun}{ixCur}.pos(1) ;%+ comparisonStructC{idxRun}{cGenAll{bIdx}.comparisonStruct{barNr}.barid}.lengthMatch;
    
    posDifMatch  = posStart(barNr);
    if comparisonStructC{idxRun}{cGenAll{bIdx}.comparisonStruct{barNr}.barid}.or(1)==2
        posDifMatch  = length(cGenAll{1}.rawBarcode)-posDifMatch-length(barcodeGen{cGenAll{bIdx}.comparisonStruct{barNr}.barid}.rawBarcode); % todo: add rescale for more accurate position
    end
    
    posD(barNr) = posTHR - (posDifExp + posDifMatch); % difference might be due to different length re-scaling factor w.r.t. theory
end

posD

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
idxPair = 3; % 1 is the pair with best score
 


import Core.plot_match_pcc;
[sortedVals, sortedIds,pscores] = sorted_scores(oScons{idFold});


[curSink,curSource] = ind2sub(size(oScons{idFold}),sortedIds(idxPair));

import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(cGen, oScons{idFold},curSink,curSource);


