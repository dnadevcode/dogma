%% Pipeline for the assembly of circular bacterial genome
% for synth data: barcoding_islands_synth_da
close all;
clear;

addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/lldev/'))
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/bargroupingprototype/'))
% addpath(genpath('/proj/snic2022-5-384/users/x_albdv/reps/hca/'))

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sF = 0.95:0.01:1.05; % length re-scaling factor
minOverlap = 300; % minimum overlap

%% Step 1: generate kymograph/mp data via lldev's genome_pipeline. Either direct kymos (from old dbm) or movies
userDir = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/ecoli_2_new/12_20_2/';

import OldDBM.General.SettingsWrapper;
defaultSettingsFilepath = SettingsWrapper.get_default_newDBM_ini_filepath();
if not(exist(defaultSettingsFilepath, 'file'))
    defaultSettingsFilepath = '';
end

dbmOSW = SettingsWrapper.import_dbm_settings_from_ini(defaultSettingsFilepath);
dbmOSW.DBMSettingsstruct.minLen = 150;
dbmOSW.DBMSettingsstruct.minOverlap = 150;

import DBM4.GenomAs.run_genome_assembly_pipeline;
[barcodeGen,barGenMerged,kymoStructs] = run_genome_assembly_pipeline(userDir, dbmOSW);

%% Step 2: Calculate overlaps (default bitmasks, perhaps succeptible to false positives)
tic
bars = barGenMerged(cellfun(@(x) sum(x.rawBitmask),barGenMerged)>minOverlap); % only for those larger than min overlap
[oS] = calc_overlap_mp(bars,sF, minOverlap,timestamp);
toc
save(fullfile(userDir,'_sessiondata','mp_data.mat'),'minOverlap','oS','sF')
barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),bars,'un',false);...
    cellfun(@(x) x.rawBitmask,bars,'un',false)]',{'rawBarcode','rawBitmask'},2);

%% Potential additions:
%%Similarity test v1: calculate std for N best match. This will be used to filter out
% % bars locally matching to several.
% import Nullmodel.local_sim_test;
% [curStd,pairs,scrs] = local_sim_test(pscores, barcodeGen, overlapStruct,minOverlap,150,3000,timestamp,overlapStructMP);
% %%%%
%% Check for barcodes which are periodic: run self similarity
% % import Core.find_periodic_barcodes;
% % [valsPassingThresh,pccValsMax,mp1,mpI1] = find_periodic_barcodes(barStruct,timestamp)

% barStruct = barStruct(find(valsPassingThresh));


%% Step 3: min length
fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/DA32087.fasta'};
% fastaFile = {'/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/018_final_polish.fasta'};
addpath(genpath('/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/'))


minLen = [200:50:2000];
% [MP,mpMax,theoryStructRev,MPI] =
% bargrouping_minimum_length([],110,0.2,1,30,[200], 5*10^6); % if thry not
% know
[MP,mpMax,theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,110,0.25,1,30,minLen, 0);


%% Step 4: good barcodes pairs
% find the mp thresh for local similarity (theory against itself)
% mpThresh = mpMax(minLen==minOverlap);
thresCC = 0.9;
import Core.filter_good_scores
[sortedVals,sortedIds,pscores] = filter_good_scores(oS,mpMax,minLen,thresCC);

% plot the graph
import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(size(oS),sortedIds,length(pscores),barGenMerged,zeros(1,length(pscores)),100,oS);

%% Step 5: greate barcode island
% pass already sorted values and indices to create_barset
import Core.create_barset_os; 
[barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = ...
    create_barset_os(sortedVals,sortedIds, barStruct,oS);

nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands)); % keep only non-empty
barcodeIslandsData = barcodeIslandsData(nonEmpty);
barcodeIslands = barcodeIslands(nonEmpty);
barislandLength = cellfun(@(x) max(x(:,2))-min(x(:,1)),barcodeIslandsData)


% plot group quick
sourceGroup = 2;
plot_temp(barcodeIslandsData{sourceGroup}, barcodeIslands{sourceGroup},barStruct )


% assembly w/o reference / bars are shuffled..
nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
import Plot.islandsPosStruct;
posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
import Plot.plot_best_pos;
fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6


[outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(bars(cellfun(@(x) x.barid,posStruct))',posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp,1);


import Plot.plot_island;
plot_island(posStruct,islandElts, bars(cellfun(@(x) x.barid,posStruct))',1);

%% Step 6: compare to theory
%% comp to thry 
fastaFile =[]; % todo:test

fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/DA32087.fasta'};
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};


nmbp = 0.236;
nmpx = 110;
psffac = 1;
nmbp = nmbp*psffac;
nmpx = nmpx*psffac;

import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmpx);

sets.comparisonMethod = 'mass_pcc';

% sF = 0.9:0.01:1.1;
[comparisonStruct] = compare_to_t(barGenMerged,theoryStruct,sF,sets)


bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

super_quick_plot(46,barGenMerged,comparisonStruct,theoryStruct)


[compStr] = compare_to_t_mp(bars,theoryStruct,sF,300)

bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], compStr, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6


%% reference based assembly: TODO: possible improvements?: for this: compare against thr
[outConsensus, coverage, pval] = gen_reference_based_assembly(barGenMerged,comparisonStruct,theoryStruct,timestamp);
% [outConsensus, coverage] = gen_assembly(barcodeGen,comparisonStruct,cumsum(repmat(10000,1,1)),timestamp);

% assembly w/o reference / bars are shuffled..
% [outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(barcodeGen(cellfun(@(x) x.barid,posStruct)),posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp);

%% Mask 'bad' regions
% 
% barcodes = arrayfun(@(x) barStruct(x).rawBarcode,1:length(barStruct),'un',false);
% bitmasks = cellfun(@(x, y) y & abs(x - nanmean(x(y))) <= 2.5 * nanstd(x(y)),...
%     arrayfun(@(x) barStruct(x).rawBarcode,1:length(barStruct),'un',false),...
%     arrayfun(@(x) barStruct(x).rawBitmask,1:length(barStruct),'un',false), 'un', 0);

%% visualize true output
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], synthStr, [], [], [],theoryStruct{1}.length,1);%1/bpPx*10^6

%% visualize best aligment
% idxPair= 1;
% import Core.plot_match_pcc;
% [sortedVals, sortedIds,pscores] = sorted_scores(synthStr2);
% [curSink,curSource] = ind2sub(size(synthStr2),sortedIds(idxPair));
% [f] = plot_match_pcc(barStruct, synthStr2,curSink,curSource,barStruct);

%% create barcoding islands
% import Core.create_barset;
% [barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = create_barset(pscores,sum(~isnan(pscores(:))), barcodeGen, synthStr2,synthStr);

%% plot barcode islands
% nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
% import Plot.islandsPosStruct;
% posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], synthStr(cell2mat(barcodeIslands)), [], [], [],10000,1);%1/bpPx*10^6



%% plot bargrouping generated island:
% pthresh = 0.00001;
% import Plot.plot_island;
% import Core.simple_bargroup;
% plot_island(posStruct,islandElts, barcodeGen(cellfun(@(x) x.barid,posStruct)),1);

% %% calculate scores:
% sF = 1;
% tic
% import Core.calc_overlap_pcc_sort_m; % this can be used if barcodes are relaxed
% [overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap,1);
% toc
%     
% res.overlapStruct = overlapStruct;
% 
% % timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% mkdir(strcat('output',timestamp));
% 
% to check local similarity/remove false positives. Skip for now for speed
% [oS] = calc_overlap_mp(barcodeGen',sF, minOverlap-50,timestamp)



%% Pairwise compare the barcodeIslands
% assembly w/o reference / bars are shuffled..
import Plot.islandsPosStruct;
posStruct = islandsPosStruct(barcodeIslandsData, barcodeIslands);

[outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(bars(cellfun(@(x) x.barid,posStruct))',posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData))),timestamp,0);

barIslandsGen = cell(1,length(consensus));
for i=1:length(consensus)
    barIslandsGen{i}.rawBarcode = consensus{i};
    barIslandsGen{i}.rawBitmask = ~isnan(consensus{i});
end
[oSIslands] = calc_overlap_mp(barIslandsGen,sF, 300,timestamp);


barStructIs = cell2struct([cellfun(@(x) double(x.rawBarcode),barIslandsGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barIslandsGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

import Core.plot_match_pcc;
[sortedValsIs, sortedIdsIs,pscoresIs] = sorted_scores(oSIslands);

lastTrueOverlap = find(sortedVals<mpMax(5),1,'first');
%%
idxPair =   1;

[curSink,curSource] = ind2sub(size(oSIslands),sortedIdsIs(idxPair));

import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(barStructIs, oSIslands,curSink,curSource);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE PLOTS TO VISUALLY CHECK IF ALG IS WORKING
% 
%% StepS1: Plot all kymos 
figure,tiledlayout(ceil(sqrt(length(kymoStructs))),ceil(length(kymoStructs)/sqrt(length(kymoStructs))),'TileSpacing','none','Padding','none')
for i=1:length(kymoStructs)
    nexttile 
        imshowpair(imresize(kymoStructs{i}.alignedMask,[200 500]),imresize(kymoStructs{i}.alignedKymo,[200 500]), 'ColorChannels','red-cyan'  )
%     imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  )
    title(num2str(i));
end

%% Plot kymo length change


%% Plot best overlap
barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),bars,'un',false);...
    cellfun(@(x) x.rawBitmask,bars,'un',false)]',{'rawBarcode','rawBitmask'},2);
% 
idxPair =   10;
[sortedVals, sortedIds,pscores,fullScores] = sorted_scores(oS);
[curSink,curSource] = ind2sub(size(oS),sortedIds(idxPair));
import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, oS,curSink,curSource);


%% Plot same but theory positions
bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct([curSink curSource]), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6
fig1 = plot_best_pos([], compStr([curSink curSource]), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

super_quick_plot(curSink,barGenMerged,comparisonStruct,theoryStruct)
super_quick_plot(curSource,barGenMerged,comparisonStruct,theoryStruct)


% super_quick_plot(curSink,barGenMerged,compStr,theoryStruct) % works for
% non-mp only
% super_quick_plot(curSource,barGenMerged,compStr,theoryStruct)

%% mol length diff. put into lldev fun
[sortedVals, sortedIds,pscores,fullScores] = sorted_scores(oS);

import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(pscores,1000,barGenMerged,zeros(1,length(pscores)),100,oS);

%%

figure,tiledlayout(ceil(sqrt(length(kymoStructs))),ceil(length(kymoStructs)/sqrt(length(kymoStructs))),'TileSpacing','none','Padding','none')
for i=1:length(kymoStructs)
    nexttile;        imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  );    title(num2str(i));
% imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  )
end