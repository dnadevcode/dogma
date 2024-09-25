% Required steps to generate results for synthetic data

% Figures:
%   Fig3

% Load pre-calculated data to save time.
% This is data generate via  figs_validation_on_synthetic_data:
% load('synth_2023-11-28_09_43_54bml.mat');
% load('synth_2023-11-28_09_43_54resRun.mat');

% load('synth_2024-06-13_22_23_19resRun.mat') % random sequence
% load('synth_2024-06-13_22_23_19resRun.mat') % random sequence
% load('synth_2024-06-13_22_23_19resRun.mat') % random sequence
% load('synth_2024-06-14_09_07_42resRun.mat')

load('synth_2024-09-05_09_03_37resRun.mat'); % nc_000913.3

%% Fig3: Plot with graphs of barcode islands
sets.minOverlap = 300; % minimum overlap
snrv = setsGen.snvr;

final_figure_3(sets, resRun, bG,snrv)


%     %% FigS10
% import Validation.bml_validation;
% [tpr, fpr] = bml_validation(synthStr,bG,theoryStruct,N, pxpsf, mL, stdL, ...
%     snr, stdE, isC, tL, minL );


%% Fig4: Plot for validation:
final_figure_4(sets, resRun, bG,synthStr,theoryStruct)

% Fig5
final_figure_5(resRun,bG,theoryStruct)

    %% make sure that above also calcc for this
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

ixtest = 5;

nB = 200;

sets.minOverlap = 300;
oS = resRun{ixtest}.oS;
barcodeGen = bG{ixtest}';

oS = oS(1:nB,1:nB);
barcodeGen = barcodeGen(1:nB);

% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, 0.08,0.01); %localCCStruct


sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.04;
sets.pxDifSetting = 50;

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],1)
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

    import Core.create_graph_from_barsetgen;
 [Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen,[],1);




tS = theoryStruct{ixtest};
bT =  theoryStruct{ixtest};

synthS = synthStr{ixtest};

    meanSF = mean(barcodeIslandsData{1}(:,end));

% approxNMBP = meanSF*nmbp

% validation against theory
import Validation.barislands_vs_theory_validation_synth

% barislands_vs_theory_validation_synth

%% TODO: rescale barcodeIslandsData{1} to a specific length re-scaling factor

% this function updates data vector. Still need to get barcodes for this
import Core.update_sf_barset;

idx = 3;

[newData] = update_sf_barset(barcodeIslandsData{idx}, 1,1);
% [newData2] = update_sf_barset(newData, 1,2) % reverse, to check that

% [newData2] = update_sf_barset(newData, 1/1.1) % reverse, to check that
% it's correct


bars = barcodeGen(barIslands{idx}); 

barx = barIslands{idx};

%todo: unique sf for each bar

sets.comparisonMethod = 'mass_pcc';


barsC = bars;
% convert bars
for i=1:length(barsC)
    lenBarTested = length( barsC{i}.rawBarcode );
     barsC{i}.rawBarcode = interp1(barsC{i}.rawBarcode, linspace(1,lenBarTested,lenBarTested*newData(i,4)));
     barsC{i}.rawBitmask =  barsC{i}.rawBitmask (round(linspace(1,lenBarTested,lenBarTested*newData(i,4))));

     if newData(i,3)==2
         barsC{i}.rawBarcode = fliplr(barsC{i}.rawBarcode);
         barsC{i}.rawBitmask = fliplr(barsC{i}.rawBitmask);
     end         
end

% compare to theory: either pcc or mp
[comparisonStruct,rezMax,~] = compare_to_t(bars,tS,sF,sets);

% can find pos-dif histogram for individual barcodes
posTrue = arrayfun(@(x) synthStr{5}{x}.pos, barx);
posFound = cellfun(@(x) x.pos(1), comparisonStruct);

figure,plot(posTrue-posFound)

% so individual ones are mapped very nicely. 

%% Now, let's fix one of the individual barcodes and recover the found
% positions from it
ix = 55 ;

bbS = comparisonStruct{ix}.bestBarStretch;  % best bar stretch (bbS)
bbO = (comparisonStruct{ix}.or(1)~=barcodeIslandsData{idx}(ix,3))+1;

% idx=1

% bars = barcodeGen(barIslands{idx}); 
% 
% barx = barIslands{idx};


[bbSAll] = update_sf_barset(barcodeIslandsData{idx}, bbS, bbO);


curBestPosFound = comparisonStruct{ix}.pos(1);

bestPosFound = bbSAll(:,1)+curBestPosFound-bbSAll(ix,1);

posTrue = arrayfun(@(x) synthStr{5}{x}.pos, barx);

pDif = posTrue-bestPosFound';
figure,plot(pDif)
figure,histogram(pDif)


%% full compare all
import Core.update_sf_barset;

idx = 3;

bars = barcodeGen(barIslands{idx}); 

barx = barIslands{idx};

sF = 0.9:0.01:1.1;

for j=1:length(sF)
    [newData] = update_sf_barset(barcodeIslandsData{idx}, sF(j),1);
    barsC = bars;
    % convert bars
    for i=1:length(barsC)
        lenBarTested = length( barsC{i}.rawBarcode );
         barsC{i}.rawBarcode = interp1(barsC{i}.rawBarcode, linspace(1,lenBarTested,lenBarTested*newData(i,4)));
         barsC{i}.rawBitmask =  barsC{i}.rawBitmask (round(linspace(1,lenBarTested,lenBarTested*newData(i,4))));
    
         if newData(i,3)==2
             barsC{i}.rawBarcode = fliplr(barsC{i}.rawBarcode);
             barsC{i}.rawBitmask = fliplr(barsC{i}.rawBitmask);
         end         
    end
    
    distRes = [];
    for i=1:length(barsC)
        [maxcoef, pos, or, secondPos, lenM,distRes{i}] = masked_MASS_PCC(tS{1}.rawBarcode, barsC{i}.rawBarcode, logical(barsC{i}.rawBitmask),true(1,length(tS{1}.rawBarcode)),2^15,0,0);
    end

end



%%
import Core.create_consensus
[consensusBar] = create_consensus(barcodeIslandsData{idx}, bars);
  



% [comparisonStruct,rezMax,~] = compare_to_t(bars,tS,sF,sets);
