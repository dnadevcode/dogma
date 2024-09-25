% Schematics of de-novo assembly problem using DNA barcodes.
rng('default')

snr = 5; 
import Core.load_synth_data;
% [bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(30,2.72,850,100,snr,0.05,1,2000,150);
[bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(50,2.72,850,100,snr,0.05,1,2000,150);


% we have pre-generated, so just load pre-generated data
% load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/synth_2023-11-28_09_43_54bml.mat');
% load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/synth_2023-11-28_09_43_54resRun.mat');

sets.minOverlap = 300;
timestamp ='';
idxRun = 1;
oS = synthStr2{idxRun};%resRun{idxRun}.oS;
barcodeGen = bG{idxRun}';

% both on local and global
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver,pvalCombined,  sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,nan, 0.08,0.01); %localCCStruct

sortedVals = sortedVals(1:length(localScore));
% sortedVals(100:end) = nan; % keep less scores to show two islands
% sortedVals(2000:end) = nan; % keep less scores to show two islands
sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

sets.scDiffSetting = 0.02;
sets.pxDifSetting = 20; %20

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] = ...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],1);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

% import Core.create_graph_from_barsetgen;
% [Gtemp] = create_graph_from_barsetgen(barsetGen,barcodeGen);

% % get best indices
% consSize = cellfun(@(x) size(x,1),outConsensus);
% 
% [sV,sI] = sort(consSize,'descend');

% this function updates data vector. Still need to get barcodes for this
import Core.update_sf_barset;
% barIslands =barIslands(1);
ix = 1;
 
pDif = cell(1,length(barIslands));
sorti =  cell(1,length(barIslands));
val =  cell(1,length(barIslands));
for idx = 1:length(barIslands)

    
    bars = barcodeGen(barIslands{idx}); 
    barx = barIslands{idx};
    
    synCur = synthStr{1}(barx);


    bbS = synCur{ix}.bestBarStretch;  % best bar stretch (bbS)
    bbO = (synCur{ix}.or(1)~=barcodeIslandsData{idx}(ix,3))+1;


    [bbSAll] = update_sf_barset(barcodeIslandsData{idx}, bbS/barcodeIslandsData{idx}(ix,4), bbO);

    % need to also update positions for all the synCur barcodes based on
    import Core.synth_to_table;
    [tableS] = synth_to_table(synCur);

    % update synth table. Want to have the same
    [tableSUpd] = update_sf_barset(tableS, tableS(ix,4)/bbS, 1);


    bbSAll(:,1:2) = bbSAll(:,1:2)-bbSAll(ix,1)+tableSUpd(ix,1); %
   
    bestPosFound = bbSAll(:,1);
    
    posTrue = tableSUpd(:,1);
    %arrayfun(@(x) synthStr{1}{x}.pos, barx);
    
%     pDif{idx} = min([sqrt((posTrue'-bestPosFound).^2) sqrt((posTrue'+theoryStruct{1}{1}.length-bestPosFound).^2)...
%         sqrt((posTrue'-theoryStruct{1}{1}.length-bestPosFound).^2)]');


    pDif{idx} = posTrue-bestPosFound;

     pDif{idx}( pDif{idx}<-theoryStruct{1}{1}.length) =   pDif{idx}( pDif{idx}<-theoryStruct{1}{1}.length) +theoryStruct{1}{1}.length;
     pDif{idx}( pDif{idx}>theoryStruct{1}{1}.length) =   pDif{idx}( pDif{idx}>theoryStruct{1}{1}.length) -theoryStruct{1}{1}.length;

    [val{idx},sorti{idx}] = sort(posTrue);

end

figure,tiledlayout(length(pDif),1)
nm = {'(D)','(E)','(F)'}
for idx=1:length(pDif);
nexttile([1 1])
hold on
plot(val{idx},     pDif{idx}(sorti{idx}),'blackx')
    xlabel('Position (px)')
    ylabel('pd')

    title([nm{idx},'PD barcode island (',num2str(idx),')'], 'Interpreter','latex')

end
% Maybe 