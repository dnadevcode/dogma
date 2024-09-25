function [] = length_validation()
%% length validation on synthetic data. For figure S11
idxRun = 1;
rng("default");

snvratio = 1;
N = 10000;
import Core.load_synth_data;
[bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(50,2.72,850,100,snvratio,0.05,1,N,150,0,[]);
barcodeGen = bG{idxRun};


% minL. Alternative is match score
% minLen = 100:50:2000;
nmPerPx = 110;
nmbp = 0.2; % length re-scale factor

% testing minOverlapWindow
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% sF = 0.85:0.025:1.15; % length re-scaling factor
sF = 0.9:0.025:1.1; % length re-scaling factor

minOverlapList = 100:50:600; % minimum overlap
minOverlapContig = 600;

alphanu = 0.09;
alphaN = 0.42;
betanu = 0.085;
betaN = 0.004;

pthresh = 0.05;

plotfigs = 0;

% import Core.calc_pos_dif;

%%
goodBars = zeros(1,length(minOverlapList));
posShiftBad = zeros(1,length(minOverlapList));
for idxO = 1:length(minOverlapList)
    minOverlap = minOverlapList(idxO);

    lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    tic
    % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
    [oS] = calc_overlap_mp(barcodeGen(lengths>=minOverlap)',sF, minOverlap,timestamp);
    toc

    % both on local and global
    [sortedVals, sortedIdsAll,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS,minOverlap, alphanu,pthresh,[],alphaN,betanu,betaN); %localCCStruct
    goodBars(idxO) = sum(~isnan(sortedVals));

    sortedIdsGood = sortedIdsAll(~isnan(sortedVals));
    sortedVals = sortedVals(~isnan(sortedVals));
    % calculate true positives:
    if ~isempty(sortedIdsGood)

        import Validation.calc_distance_scaled; % need to properly scale things when calculating distances for validation
        posShift = calc_distance_scaled(oS,sortedIdsGood, synthStr{1}(lengths>=minOverlap),theoryStruct{1}{1}.length, length(sortedVals));
    
        posShiftBad(idxO) = sum(posShift > 30 );
    end



end

%%
numOverlapsAbove = arrayfun(@(x) sum([synthStr2{1}.overlaplen]>x),minOverlapList);
% numOverlapsAbove = arrayfun(@(x) sum([synthStr2{1}.overlaplen]>x),50:50:1500);

f = figure
% f = figure;tiledlayout(2,1)
% nexttile
% plot(minOverlapList,goodBars./numOverlapsAbove(1))
hold on
% plot(minOverlapList,goodBars./numOverlapsAbove)
plot(minOverlapList,goodBars)

plot([300 300],[0 max(goodBars)],'red-')
xlabel('Local alignment length','Interpreter','latex');
ylabel('Number of significant overlaps','Interpreter','latex')
legend({'Number of significant overlaps','Default overlap length w'},'location','southoutside','Interpreter','latex')
print('FIGS/FigS5.eps','-depsc','-r500');

%TODO: maybe also plot number of false positives

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
print('FIGS/FigS5.eps','-depsc','-r300');
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



end

