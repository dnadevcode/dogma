%% Pipeline for the assembly of circular bacterial genome
% barcoding_islands_synth_bacterial for the synthetic case

nmPerPx = 110; % nm/px (can be extracted from image info file)
nmbp = 0.2; % nm/bp (can be extracted via dna_barcode_matchmaker lambda pipeline
psffac = 1; % scaling factor (in case PSF needs to be something else than 300nm)
numWorkers = 30; % num workers, in one node there are 30
minLen = 200; % should be same as minOverlap 
% lenSeq = 5*10^6; % approximate length of genome
lenSeq =5174317;

% clear and close all
close all;
clear all;

%% Loading data and calculating local alignment
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

% load barcode data / output from dna_barcode_matchmaker
data = {'/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/1219/_sessiondata/session_data.mat',...
    '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/12_20_1/_sessiondata/session_data.mat',...
    '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/12_20_1/_sessiondata/session_data.mat'};
%

fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/DA32087.fasta'};
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

idx = 1;

try % loading the data
    load(data{idx}); % loads barcodeGen and barcodeMerged
    load(strrep(data{idx},'session_data.mat','mp_data.mat')); % loads oS, sF
catch
    % otherwise need to run dna_barcode_matchmaker on this data
    %     dna_barcode_matchmaker
end
%

%% barcode similarity
%% If theory bar not know, run this for random theory (of expected length) against itself
minLen =[200:50:1000];
[MP,mpMax] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,psffac,numWorkers,minLen, lenSeq)
[MPS,mpMaxS] = bargrouping_minimum_length([],nmPerPx,nmbp,psffac,numWorkers,minLen, lenSeq)


% figure, errorbar(minLen,mean(cell2mat(mpMax')),std(cell2mat(mpMax')))
figure
plot(minLen,mpMaxS)
hold on
plot(minLen,mpMax)
xlabel('Overlap length')
ylabel('Max overlap PCC')
legend({'Simulated','E-coli'})


% find the mp thresh for local similarity (theory against itself)
mpThresh = mpMax(minLen==minOverlap);

%
minLen = 200; % minimum length
numF = 10; % number of frames 
sF = 1;%0.95:0.01:1.05;
minOverlap = 250;
scorethresh = 0.7; % instead use Stouffer like for HC

    
matFile = 1;
if matFile
    load(matFileLoc);
else
    barcodeGen = barcodeGen1;
end


barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

%% Check for barcodes which are periodic: run self similarity

% % import Core.find_periodic_barcodes;
% % [valsPassingThresh,pccValsMax,mp1,mpI1] = find_periodic_barcodes(barStruct,timestamp)

% barStruct = barStruct(find(valsPassingThresh));

%% check and merge similar neighbours (in case this was not identified via DBM)

%% Mask 'bad' regions

barcodes = arrayfun(@(x) barStruct(x).rawBarcode,1:length(barStruct),'un',false);
bitmasks = cellfun(@(x, y) y & abs(x - nanmean(x(y))) <= 2.5 * nanstd(x(y)),...
    arrayfun(@(x) barStruct(x).rawBarcode,1:length(barStruct),'un',false),...
    arrayfun(@(x) barStruct(x).rawBitmask,1:length(barStruct),'un',false), 'un', 0);

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


%% comp to thry 
fastaFile =[]; % todo:test

fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};


% nmbp = 0.216;
nmpx = 110;
psffac = 1;
nmbp = nmbp*psffac;
nmpx = nmpx*psffac;

import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmpx);

sets.comparisonMethod = 'mass_pcc';

sF = 0.9:0.01:1.1;
[comparisonStruct] = compare_to_t(barcodeGen,theoryStruct,sF,sets)


bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

% MP:
sF = 0.95:0.01:1.05;
[compStr] = compare_to_t_mp(bars,theoryStruct,0.95,200)

bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], compStr, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6


%% reference based assembly: TODO: possible improvements?: for this: compare against thr
[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,comparisonStruct,theoryStruct,timestamp);
[outConsensus, coverage] = gen_assembly(barcodeGen,comparisonStruct,cumsum(repmat(10000,1,1)),timestamp);

% assembly w/o reference / bars are shuffled..
[outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(barcodeGen(cellfun(@(x) x.barid,posStruct)),posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp);

%% plot bargrouping generated island:
% pthresh = 0.00001;
% import Plot.plot_island;
% import Core.simple_bargroup;
% plot_island(posStruct,islandElts, barcodeGen(cellfun(@(x) x.barid,posStruct)),1);

%% calculate scores:
sF = 1;
tic
import Core.calc_overlap_pcc_sort_m; % this can be used if barcodes are relaxed
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap,1);
toc
    
res.overlapStruct = overlapStruct;

% timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
mkdir(strcat('output',timestamp));

% to check local similarity/remove false positives. Skip for now for speed
% [oS] = calc_overlap_mp(barcodeGen',sF, minOverlap-50,timestamp)

%% now barset - calculate p-values.  


[sortedVals, sortedIds,pscores] = sorted_scores(oS);
lastTrueOverlap = find(sortedVals<mpThresh,1,'first');


% import Zeromodel.pval_for_struct;
% pscores = pval_for_struct(overlapStruct,scorethresh,0.05);
% 
% 
% psc = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1
% filtM = triu(psc);
% filtM(filtM==0) = nan;

% lastTrueOverlap = 100;
% pass already sorted values and indices to create_barset
import Core.create_barset_os; %circular case?
[barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = ...
    create_barset_os(sortedVals(1:lastTrueOverlap-1),sortedIds(1:lastTrueOverlap-1), barStruct,oS);

sourceGroup = 11;
plot_temp(barcodeIslandsData{sourceGroup}, barcodeIslands{sourceGroup},barStruct )

nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));


barcodeIslandsData = barcodeIslandsData(nonEmpty);
barcodeIslands = barcodeIslands(nonEmpty);

barislandLength = cellfun(@(x) max(x(:,2))-min(x(:,1)),barcodeIslandsData)

% assembly w/o reference / bars are shuffled..
nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));
import Plot.islandsPosStruct;
posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
import Plot.plot_best_pos;
fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6


[outConsensus, coverage, consensus,islandElts, islandPx] = gen_assembly(bars(cellfun(@(x) x.barid,posStruct))',posStruct,cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),timestamp,1);


import Plot.plot_island;
plot_island(posStruct,islandElts, bars(cellfun(@(x) x.barid,posStruct))',1);

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
idxPair =   7;

[curSink,curSource] = ind2sub(size(oSIslands),sortedIdsIs(idxPair));

import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(barStructIs, oSIslands,curSink,curSource);

