function [oS,barcodeGen, barsetGen, outConsensus, coverage, consensus, ...
    islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands,sortedVals, sortedIds] = prep_data(dataFold,idx)

    % prep_data - needs to be updated with correct parameters!
    
if nargin < 2
    idx = 1;
end
switch dataFold

    case 'dEC'
        % load data. Here it changes to folder with external data. For
        % reproducibility, the '.mat' file with the overlaps can be downloaded
        cd('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/')
        load('3all_2023-11-24_12_38_55.mat')
        %     load('3_2023-11-22_11_34_44resRun.mat')
        %     load('3_individualdays.mat')

        %% Parameters
        alphaNu = 0.075; % 0.08 default for the paper


    case 'dECind'
        cd('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/')
        %         load('3all_2023-11-24_12_38_55.mat')
        load('3_2023-11-22_11_34_44resRun.mat')
        oS = resRun{idx}.oS;
        barcodeGen = bG{idx};
        %     load('3_individualdays.mat')

        %% Parameters
        alphaNu = 0.075; % 0.08 default for the paper
    otherwise

end

sets.minOverlap = 300;
sets.scDiffSetting = 0.04;
sets.pxDifSetting = 50;
showFig = 0; % whether to show output figure
pthresh =  0.01;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'); % timestamp for results

%%
% both on local and global. Only uses the oS structure as input +
% parameters from table 2
[sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
    calculate_sorted_pvals(oS,sets.minOverlap, alphaNu,pthresh); %localCCStruct

sortedValsGood = sortedVals(~isnan(sortedVals));
sortedIdsGood  = sortedIds(~isnan(sortedVals));
localScore = localScore(~isnan(sortedVals));

% sortedValsGood(4190:end) = nan; % change alphaNu instead of this

% Main barcode islands script:

import Core.barcode_island_output;
[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],showFig);
%     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

% nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty
accepted = cellfun(@(x) x.accepted, barsetGen.currentData); % all the

import Core.create_graph_from_barsetgen;
[Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen,[],showFig);

end

