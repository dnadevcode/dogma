% version of chrom_assembly_script_vs_ref_synthetic
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% STEP1: gen synthetic data

%% gen synthetic
import Nullmodel.gen_random_fragments;

% parameters
TOTAL_RAND_LENGTH = 10000; % total length
PSF_WIDTH_PIXELS = 2.7; % psf
NUM_RAND_FRAGMENTS = 50; % number of random fragments
MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
FRAGMENT_LENGTH_STD = 100; %std of length of fragment
ADDED_NOISE_MEAN = 0.4; % additive noise mean
ADDED_NOISE_STD = .1; % additive noise std
FRAGMENT_STRETCH_STD = 0.02; % fragment length-rescale factor std 0.02
IS_CIRCULAR = 1; % whether circular barcode
    
% generate barcodes
[barcodeGen,refBarcode,origPos, origFlip,origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD,IS_CIRCULAR);

% barSynth = cell2struct([barcodes'; bitmasks']',{'rawBarcode','rawBitmask'},2);
barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);

% lengths of barcode masks
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);


% plot the fragment positions. Requires hca function
synthStr = [];
for i=1:length(barcodeGen)
    synthStr{i}.idx = 1;
    synthStr{i}.pos = origPos(i);
    synthStr{i}.lengthMatch = round(length(barcodeGen{i}.rawBarcode)*origStr(i));
    synthStr{i}.or = origFlip(i);
    synthStr{i}.rf = origStr(i);
end
bpPx = 361;
import CBT.Hca.UI.Helper.plot_best_pos;
plot_best_pos([], synthStr, [], [], [],TOTAL_RAND_LENGTH,1/bpPx*10^6);

%

import Nullmodel.plot_random_fragments
[outConsensus] = plot_random_fragments(barcodeGen,1:length(barcodeGen),refBarcode,origPos,1./origStr,origFlip'+1);

sF = [0.9:0.002:1.1];



[a,msg] = mkdir(strcat('output',timestamp));
i=1;theoryStruct=[];
theoryStruct{i}.filename = fullfile(strcat('output',timestamp),'theory_barcode.txt');
fileID = fopen(theoryStruct{i}.filename,'w');
fprintf(fileID,strcat([' %2.' num2str(14) 'f ']), refBarcode);
fclose(fileID);
theoryStruct{i}.meanBpExt_nm = 0.3;
theoryStruct{i}.psfSigmaWidth_nm = 300;
theoryStruct{i}.length = length(refBarcode);
theoryStruct{i}.isLinearTF = 0;
theoryStruct{i}.name = 'Synthetic theory';

%% STEP2: MP vs theoretical barcode
% re-scale theoryGen)
sF = 0.95:0.01:1.05;
MIN_OVERLAP_PIXELS = 300;
numWorkers = 32;

% sF = 1;

minLen = 300; %150kb? or less? depends on application // if GUI, user selects thi

ix = 1;
import Core.compare_mp_all;
[cs{ix}.mpI1,cs{ix}.mp1,cs{ix}.maxMP,cs{ix}.stridx,cs{ix}.compStr] = compare_mp_all(theoryStruct,barcodeGen',minLen,ix, timestamp,sF,MIN_OVERLAP_PIXELS,numWorkers);

lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));
% lengthBorders = lengthBorders(ix);
% fig1 = figure;
bpPx = 361;
import CBT.Hca.UI.Helper.plot_best_pos;
plot_best_pos([], cs{1}.compStr, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

tryePos = origPos;

barSt = cellfun(@(x) find(x.rawBitmask==1,1,'first'),barcodeGen);
detPos = cellfun(@(x) x.pB-x.pA+1,cs{ix}.compStr)-barSt'+1; % check more precicely.
detPos(detPos<0)=detPos(detPos<0)+theoryStruct{1}.length;

sumSq = sqrt((origPos-detPos).^2);
figure,plot(sumSq)
% figure,plot(origPos-detPos)

figure,plot(cs{ix}.compStr{16}.allposB-cs{ix}.compStr{16}.allposA)
% 
% detPos = cellfun(@(x) mean(x.allposB(1:100)-x.allposA(1:100))+1,cs{ix}.compStr)-barSt'+1; % check more precicely.
% detPos(detPos<0)=detPos(detPos<0)+theoryStruct{1}.length;
% sumSq = sqrt((origPos-detPos).^2);
% figure,plot(sumSq)

[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,cs{1}.compStr,theoryStruct);

%%
bG = barcodeGen';
minLen = 300; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),bG);
bG = bG(barLens>minLen);

MIN_OVERLAP_PIXELS = 200;
% NMPX = 110;
% NMPSF = 200;
% NMBP = 0.22;
sF = 0.9:0.01:1.1; %0.85:0.01:1.15;
% sF = 0.85:0.01:1.15;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);


sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP' ;
if ispc % for PC get from the initialization
    SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else %
    sets.SCAMP_LINE = '~/SCAMP/' ;

    SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end

%% save barcodex into txts

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,bG,'un',false);...
    cellfun(@(x) x.rawBitmask,bG,'un',false)]',{'rawBarcode','rawBitmask'},2);
% have to save separately..
foldSynth= 'barcoli';
% have to save separately..
 [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);


% skip one barcode
[names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);


%%
MIN_OVERLAP_PIXELS = 300;

numWorkers = 32;
%%  calculates all MP and MPI
tic
import Core.calc_overlap;
[mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,names2,namesBar);
toc

% tic
% [mpIBA,mpBA,maxMPBA] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,namesBar,names2);
% toc

for ii=1:length(baridx2)
    baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
end
% we want to create a nicer structure for all-to-all comparison and
% contains easily accesible data.
tic
import Core.mp_res_to_struct;
[overlapStruct] = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
toc
import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, overlapStruct,11,6);
PCC_OVERLAP = reshape([overlapStruct.fullscore], size(overlapStruct,1),size(overlapStruct,2))
overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2))
% lenOverlap = 

%%
%% filter scores based on best overlap
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');

% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 3;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));


PCC_OVERLAP(xId,yId)

import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, overlapStruct,xId,yId);
[res] = get_display_results(barcodeGen([xId yId]),[], comparisonStruct([xId yId]), theoryStruct, sets);
  import CBT.Hca.UI.Helper.plot_best_pos;
    plot_best_pos([], compStr([xId yId]), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

[f] = plot_match_simple(barStruct, overlapStruct,31,18);
[f] = plot_match_simple(barStruct, overlapStruct,18,31);

pthresh = 0.00001;
import Core.simple_bargroup;
[orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,xId,pthresh,1);


[res] = get_display_results(barcodeGen(bars),[], comparisonStruct(bars), theoryStruct, sets);

import Core.bargroup_map_simple;
bargroup_map_simple