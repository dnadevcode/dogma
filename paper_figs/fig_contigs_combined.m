% calculate assembly for a number of barcodes from different folders
%
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sF = 0.85:0.025:1.15; % length re-scaling factor
minOverlap = 300; % minimum overlap
minOverlapContig = 600;

plotfigs=0;

fastas = {'018_final_polish.fasta','DA32087.fasta'}; % two known fastas (if used)

% ix1= 10;
% ix2 = 12;
% barcodeGen = [bG{ix1} bG{ix2}];
idxRun = 1;
barcodeGen = horzcat(bG{:});

nmbp = 0.25; % what is the correct nmbp ratio for this experiment?

% perm = randperm(length(barcodeGen));
% barcodeGen = barcodeGen(perm(1:500));


% for idxRun = 1:length(bG)
%     barcodeGen = bG{idxRun};

lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(find(lengths>minOverlap));

nB = length(barcodeGen);
tic
% bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
[oS] = calc_overlap_mp(barcodeGen,sF, minOverlap,timestamp);
toc

save([num2str(dataSetIdx), 'all_',timestamp,'.mat'],'oS','barcodeGen','-v7.3');

save('out_ecoli_new_1.mat','oS','barcodeGen','-v7.3');

save('out_ecoli_sample_2.mat','oS','barcodeGen','-v7.3');

%     save('out_500_ecoli1.mat','oS','perm','-v7.3');



minLen = [150:50:3000]; % min to max % could take more points for more accurate.. 
[MP,mpMaxLenBasedC{idxRun},theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastas(2),110,nmbp,1,30,minLen, []);% no lenseq if based on fastafile

% alignment to predicted theory // here we can also have a combined score
% (best local) for a more accurate result
sets.comparisonMethod = 'mass_pcc';
[comparisonStructC{1},rezMax,bestBarStretch] = compare_to_t(barcodeGen,theoryStruct,0.9:0.01:1.1,sets);



%     barcodeGen = barcodeGen
bestBarStretch = reshape([oS.bestBarStretch], size(oS,1),size(oS,2));

selfS  = triu(bestBarStretch(1:nB,1:nB));
dat1 = selfS(selfS~=0);

figure,histogram(selfS(selfS~=0))
%%
selfS1  = triu(bestBarStretch(1:length(bG{1}),1:length(bG{1})));
dat1 = selfS1(selfS1~=0);
figure,histogram(selfS1(selfS1~=0))

selfS2  = triu(bestBarStretch(1:length(bG{ix1}),length(bG{ix1})+1:end));
dat2 = selfS2(selfS2~=0);
figure,histogram(selfS2(selfS2~=0))


selfS3  = triu(bestBarStretch(length(bG{ix1})+1:end,length(bG{ix1})+1:end));
dat3 = selfS3(selfS3~=0);
figure,histogram(dat3)

%%
% 
% % 
% % idxs = 1:length(bgAll);
% % barcodeGen = bgAll(barsPassThresh);idxs = barsPassThresh;
% 
% bgS = barcodeGen;
% 
% % [oS] = calc_scores({bgAll},300,0.95:0.025:1.05);
% 
% tic
% % bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
% [oS] = calc_overlap_mp(bgS,sF, minOverlap,timestamp);
% toc
% 
%% bargrouping_funs


import Core.plot_match_pcc;
[sortedValsAll, sortedIdsAll,pscores] = sorted_scores(oS);
idxThreshMP = arrayfun(@(x) find(x<=minLen,1,'first'),minOverlap);
coefDiffsMP = sortedValsAll-mpMaxLenBasedC{idxRun}(idxThreshMP);
barsPassThreshMP = find(coefDiffsMP > 0);

thresCC = mpMaxLenBasedC{idxRun}(idxThreshMP);
import Core.filter_good_scores
[sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMaxLenBasedC{idxRun},minLen,thresCC);

length(sortedVals)

%%
if plotfigs
%% start
idFold = 1;
idxPair = 1;
 


% allCoefs = cellfun(@(x) x.maxcoef(1),compStr);
% allLengths =300;
% 
% coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% 


%

[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));

% lenThry =cumsum(cellfun(@(x) x.length,theoryStruct));
import Plot.pair_evaluation_plot;
[f] = pair_evaluation_plot(barcodeGen, oS,curSink,curSource,comparisonStructC{idxRun}([curSink curSource]),sum(theoryStruct.length));
%% plot
print('FIGS/Fig1.eps','-depsc','-r300');

%%
passThreshPair = zeros(1,N);
N = 1000;
posShift = zeros(1,N);
for idxPair=1:N
    
    [curSink,curSource] = ind2sub(size(oS),sortedIdsAll(idxPair));
    if ismember(curSink,passthreshC{idxRun})&&ismember(curSource,passthreshC{idxRun})
        passThreshPair(idxPair) = 1;
    end
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
% lgn = legend({'Distance between estimated GT and pair overlap','Threshold for good scores'})
nexttile
plot(posShift(1:length(barsPassThreshMP)+50));hold on
plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
plot(find(~passThreshPair(1:length(barsPassThreshMP)+50)),posShift(find(~passThreshPair(1:length(barsPassThreshMP)+50))),'redx')
lgn = legend({'Distance between estimated GT and pair overlap','Threshold for good scores','Not accurate map against theory'});
lgn.Location ='southoutside';
print('FIGS/FigS3.eps','-depsc','-r300');
end
%%
% plot the graph
import Core.create_overlap_graph;
[finalgraph,Ggraphs,falseMaps] = create_overlap_graph(size(oS),sortedIds,length(sortedIds),barcodeGen,zeros(1,length(pscores)),100,oS);
print('FIGS/FigS4.eps','-depsc','-r300');

%


% [mpMax,finalgraph,Ggraphs,sortedVals,sortedIds,pscores] = bargrouping_funs(oS,bG,0.8)
% toc
[outConsensus, coverage, consensus,islandElts, islandPx,cGen,badData,...
    badDataInfo,barcodeIslandsData,barcodeIslands,barStruct] = ...
    bar_island_out(sortedVals,sortedIds, oS,barcodeGen,timestamp);
print('FIGS/Fig2.eps','-depsc','-r300');

[a,b] = sort(cellfun(@(x) size(x,1),barcodeIslandsData),'desc')
sourceGroup = b(1);
cGen = cGen(b(1:12));
plot_temp(barcodeIslandsData{sourceGroup}, barcodeIslands{sourceGroup},barStruct )
print('FIGS/Fig3.eps','-depsc','-r300');


%%
import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmPerPx);
barcodeGenT{1}.rawBitmask = logical(barcodeGenT{1}.rawBitmask);
% barcodeGenT{1}.rawBarcode = [barcodeGenT{1}.rawBarcode;

sets.comparisonMethod = 'mass_pcc';

%  add theory as last barcode
cGen =[cGen barcodeGenT];

sF = 0.9:0.01:1.1;
% [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(cGen,theoryStruct,sF,sets);
% [compStr,~,calcLengths] = compare_to_t_mp(cGen,theoryStruct,sF,300); % get oS as output

% calculate MP scores 
[oScons] = calc_overlap_mp(cGen,sF, minOverlapContig,timestamp);
 
import Core.plot_match_simple;
[f] = plot_match_simple(cGen, oScons, 1,length(cGen));
print('FIGS/Fig4.eps','-depsc','-r300');

[f] = plot_match_simple(cGen, oScons, 1,2);

import Core.plot_match_islands;
[f] = plot_match_islands(cGen, oScons, [1:12],length(cGen));
print('FIGS/Fig5.eps','-depsc','-r300');

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
super_quick_plot(15,barcodeGen,comparisonStructC{1},theoryStruct)


compare_consensuses()


end
%% Now validation for full


%%

import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(barcodeGen, oS,curSink,curSource);

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

%%--

N = 12;
[curSink,curSource] = ind2sub(size(oS),sortedIds(N))

[finalgraph,Ggraphs,falseMaps] = ...
    create_overlap_graph(size(oS),sortedIds(1:N),N,barcodeGen,zeros(1,length(pscores)),100,oS);
import Core.plot_match_simple;
[f] = plot_match_simple(barcodeGen, oS,100,409);
%%

[finalgraph,Ggraphs,falseMaps] = create_overlap_graph(size(oS),sortedIdsAll,1000,barcodeGen,zeros(1,length(sortedIdsAll)),100,oS);

N = 2000;
[outConsensus, coverage, consensus,islandElts, islandPx,cGen,badData,...
    badDataInfo,barcodeIslandsData,barcodeIslands,barStruct] = ...
    bar_island_out(sortedValsAll(1:N),sortedIdsAll(1:N), oS,barcodeGen,timestamp);