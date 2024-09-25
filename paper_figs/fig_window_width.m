
% testing minOverlapWindow

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sF = 0.85:0.025:1.15; % length re-scaling factor
minOverlapList = 100:50:600; % minimum overlap
minOverlapContig = 600;

plotfigs=0;


minLen = 100:50:2000;
ix1= 1;
% ix2 = 12;
barcodeGen = [bG{1} bG{2} bG{3} bG{4}];

nmPerPx = kymoStructs{ix1}{1}.nmpxidFold;
nmbp = kymoStructs{ix1}{1}.nmBpidFold;
fastaFile = fastas(fastaFiles(idxRun));
% can be calculated before for all nmbp
[MP,mpMaxLenBasedC{ix1},theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,1,numWorkers,minLen, []);% no lenseq if based on fastafile
%


%%
for idxRun = 1:length(minOverlapList)
    minOverlap = minOverlapList(idxRun);

    lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(barcodeGen(lengths>=minOverlap),sF, minOverlap,timestamp);
    toc
    import Core.plot_match_pcc;
[sortedVals, sortedIdsAll,pscores] = sorted_scores(oS);
idxThreshMP = arrayfun(@(x) find(x<=minLen,1,'first'),minOverlap);
coefDiffsMP = sortedVals-mpMaxLenBasedC{ix1}(idxThreshMP);
barsPassThreshMP = find(coefDiffsMP > 0);

thresCC = mpMaxLenBasedC{ix1}(idxThreshMP);
import Core.filter_good_scores
[sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBasedC{ix1},minLen,thresCC);

goodBars(idxRun) = length(sortedVals);
end
%%

% [f] = plot_match_simple(bgS, oS,curSink,curSource);


%%
if plotfigs
    %%
idFold = 1;
idxPair = 2;
 


% allCoefs = cellfun(@(x) x.maxcoef(1),compStr);
% allLengths =300;
% 
% coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% 


%

[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));

lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
% import Plot.pair_evaluation_plot;
% [f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,comparisonStructC{idxRun}([curSink curSource]),cumsum(cellfun(@(x) x.length,theoryStruct)));
% print('FIGS/Fig1.eps','-depsc','-r300');

import Core.plot_match_simple;
[f] = plot_match_simple(barcodeGen, oS,curSink,curSource);


%%
N = 1000;
posShift = zeros(1,N);
for idxPair=1:N
    
    [curSink,curSource] = ind2sub(size(oS),sortedIdsAll(idxPair));
    pos1 = comparisonStructC{idxRun}{curSink}.pos(1);
    pos2 = comparisonStructC{idxRun}{curSource}.pos(1);
    if comparisonStructC{idxRun}{curSource}.or(1)==1
        posDifTheory = pos1-pos2;
    else
        posDifTheory = -(pos1-pos2+comparisonStructC{idxRun}{curSink}.lengthMatch-comparisonStructC{idxRun}{curSource}.lengthMatch);
    end
    posDifExp = oS(curSink,curSource).pB - oS(curSink,curSource).pA;


    
    posShift(idxPair) = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);
end

%   passthreshC{idxRun}  these passed thresh
figure,tiledlayout(2,1)
nexttile
plot(posShift);hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})
nexttile
plot(posShift(1:100));hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})
print('FIGS/FigS3.eps','-depsc','-r300');
end
%%
% plot the graph
import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(size(oS),sortedIds,length(sortedIds),barcodeGen,zeros(1,length(pscores)),100,oS);
print('FIGS/FigS4.eps','-depsc','-r300');

%
% [mpMax,finalgraph,Ggraphs,sortedVals,sortedIds,pscores] = bargrouping_funs(oS,bG,0.8)
% toc
[outConsensus, coverage, consensus,islandElts, islandPx,cGen] = bar_island_out(sortedVals,sortedIds, oS,barcodeGen,timestamp)
print('FIGS/Fig2.eps','-depsc','-r300');

%%
import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmPerPx);
barcodeGenT{1}.rawBitmask = logical(barcodeGenT{1}.rawBitmask);

sets.comparisonMethod = 'mass_pcc';

%  add theory as last barcode
cGen =[cGen barcodeGenT];

sF = 0.9:0.01:1.1;
% [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(cGen,theoryStruct,sF,sets);
% [compStr,~,calcLengths] = compare_to_t_mp(cGen,theoryStruct,sF,300); % get oS as output

% calculate MP scores 
[oScons] = calc_overlap_mp(cGen,sF, minOverlapContig,timestamp);
 
import Core.plot_match_simple;
[f] = plot_match_simple(cGen, oScons, 8,10);
print('FIGS/Fig4.eps','-depsc','-r300');


% [f] = plot_match_simple(bgS, oS,curSink,curSource);


% 
% allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
% allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);
% 
% idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
% coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% 
% barsPassThresh = find(coefDiffs > 0.05);
% % barsPassThresh = find(ones(1,length(comparisonStruct)));

super_quick_plot(910,bgAll,comparisonStruct,theoryStruct)


compare_consensuses()


end
%% Now validation for full


%%

import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(bgS, oS,curSink,curSource);

bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct(idxs([curSink curSource])), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6


% plot on theory

% %% check position shift (ignoring orientation for now)
% N = 1000;
% posShift = zeros(1,N);
% for idxPair=1:N
%     
%     [curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
%     pos1 = comparisonStruct{idxs([curSink])}.pos(1);
%     pos2 = comparisonStruct{idxs([curSource])}.pos(1);
%     if comparisonStruct{idxs([curSource])}.or(1)==1
%         posDifTheory = pos1-pos2;
%     else
%         posDifTheory = -(pos1-pos2+comparisonStruct{idxs([curSink])}.lengthMatch-comparisonStruct{idxs([curSource])}.lengthMatch);
%     end
%     posDifExp = oS(curSink,curSource).pB - oS(curSink,curSource).pA;
% 
% 
%     
%     posShift(idxPair) = min([dist(posDifTheory,posDifExp) dist(posDifTheory+lenThry,posDifExp) dist(posDifTheory-lenThry,posDifExp) ]);
% end
% 
% figure,plot(posShift);hold on
% plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
% legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})


%% bargrouping_funs
thresCC = mpMaxLenBased(idxThreshMP);
import Core.filter_good_scores
[sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBased,minLen,thresCC);

