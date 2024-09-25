% This script runs bargrouping using overlaps


testSet = {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220221_Sample-dEC-st10_464.93bpPERpx_0.236nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220224_Sample-dEC-st11_425.53bpPERpx_0.258nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220223_Sample-dEC-st11_519.34bpPERpx_0.211nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220225_Sample-dEC-st11_553bpPERpx_0.198nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220613_Sample DA32087-4-st2_570bpPERpx_0.193nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220610_Sample DA32087-4-st1_570bpPERpx_0.192nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/ch1_test/'};

fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% 
% 
testSet = {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/211011_sample5_539bpPERpx_0.204 nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/211101_Sample 5.2_542.83bpPERpx_0.202nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/211028_Sample 5.3_ 536.22bpPERpx_0.205nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'}
fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

idx = 4;
numF = 10;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');


import Core.create_barcodegen;
[barcodeGen, kymoStructs, sets] = create_barcodegen(testSet{idx},1,numF,1,timestamp,1);

%%
% only keep barcodes longer than minimum length.
minLen = 300; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(barLens>minLen);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);


sets.comparisonMethod = 'mass_pcc';
sets.filter = 0;
nmbp = barcodeGen{1}.nmbp;
nmpx = 110; % 208?
psf = 300; %300 nm

%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)
% nmbp = 0.22;
import Thry.gen_theoretical;
[theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);

%%
%% PCC
sF = 0.8:0.01:1.2;

rezMax=[];bestBarStretch=[];bestLength=[];
for i=1:length(theoryStruct)
    tic
    import CBT.Hca.Core.Comparison.on_compare;
    [rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGen,...
        theoryStruct{i},sets.comparisonMethod, sF,[],50,[],sets.filterSettings);
    toc
end

comparisonStructAll = rezMax;
for i=1:length(comparisonStructAll)
    for j=1:length(bestBarStretch{i})
        comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
        comparisonStructAll{i}{j}.length = bestLength{i}(j);
    end
end
import CBT.Hca.Core.Comparison.combine_theory_results;
[comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);

sets.timeFramesNr = nan;
sets.displayResults=1
sets.userDefinedSeqCushion = 0;
sets.genConsensus = 0;
import CBT.Hca.UI.get_display_results;
[res] = get_display_results(barcodeGen,[], comparisonStruct, theoryStruct, sets);
% %% MP

%%
% two things to do: 

% Edge pcc: ignore circular/include circular
% 1) add zeros to second barcode to 
import SignalRegistration.masked_pcc_corr;

shortVec = barcodeGen{1}.rawBarcode;
w1  =  barcodeGen{1}.rawBitmask; % add at least length(shortVec)-1
longVec = [barcodeGen{2}.rawBarcode zeros(1,max(length(shortVec)-1,length(shortVec)-length(barcodeGen{2}.rawBarcode)))];
w2  =  [barcodeGen{2}.rawBitmask zeros(1,max(length(shortVec)-1,length(shortVec)-length(barcodeGen{2}.rawBarcode)))];
minOverlap = 300;

iy =xId
ix = yId;
barA = barcodeGen{iy}.rawBarcode;
wA = barcodeGen{iy}.rawBitmask;
barB = barcodeGen{ix}.rawBarcode;
wB = barcodeGen{ix}.rawBitmask;

extraL = max(length(barA)-1-minOverlap,length(barA)-length(barB));
% barB = zeros(1,);
% if circular: long

barB2 = [zeros(1,extraL) barB];
wB2 = [zeros(1,extraL) wB];
tic
[ xcorrs, numElts ] = masked_pcc_corr( barA,barB2,wA,wB2,minOverlap ); %todo: include division to k to reduce mem
toc

tic % slow version:
xcorrs2 = zeros(size(xcorrs));
numElts2 = zeros(size(xcorrs));
for i=1:size(xcorrs,2)
    mask = wA.*wB2(mod(i-1:i+length(wA)-2,length(wB2))+1);
    numElts2(1,i) = sum(mask);
    xcorrs2(1,i) = zscore(barA(find(mask)),1)*zscore(barB2(mod(i-2+find(mask),length(wB2))+1),1)'/sum(mask);
end
toc
% tic
% import Core.masked_pcc_edge_comp
% [xcorrs, numbits] = masked_pcc_edge_comp(shortVec, longVec, w1, w2, minOverlap);
% toc 

%% 
% delete(gcp('nocreate'))
% parpool('threads')

delete(gcp('nocreate'))
parpool('local',32)

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);


minOverlap = 300;
tic
import Core.calc_overlap_pcc;
[overlapStruct] = calc_overlap_pcc(barcodeGen, sF,minOverlap);
toc

PCC_OVERLAP = reshape([overlapStruct.score], size(overlapStruct,1),size(overlapStruct,2));
overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));
% lenOverlap = 

%%
%% filter scores based on best overlap
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');

% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 1
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

% overlapStruct(xId,yId)
PCC_OVERLAP(xId,yId)

import Core.compute_dense_zvals

% tic
% zVal = [];
% for i=1:size(overlapStruct,1)
%     i
%     for j=1:size(overlapStruct,2)
%           zVal(i,j) = compute_dense_zvals(overlapStruct(i,j).score,...
%  overlapStruct(i,j).lenA, overlapStruct(i,j).lenB, overlapStruct(i,j).h, 1, 'minimum');
%     end
% end
% toc

%         zVal(i,j) = compute_dense_zvals(overlapStruct(xId,yId).score,...
%  overlapStruct(xId,yId).lenA, overlapStruct(xId,yId).lenB, overlapStruct(xId,yId).h, 1, 'minimum');
% toc

import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,xId,yId);

% pthresh = 0.00001;
% import Core.simple_bargroup;
% [orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,yId,pthresh,1);

%% Same with MP
MIN_OVERLAP_PIXELS=300;
barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
% have to save separately..
foldSynth= 'barcoli';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);
% skip one barcode
[names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);



if ispc % for PC get from the initialization
    sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP' ;
    SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else %
%     sets.SCAMP_LINE = '/home/albyback/postdocData/test_transloc/SCAMP/';
    sets.SCAMP_LINE = '~/SCAMP/';

    SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end

delete(gcp('nocreate'))
numWorkers = 32;
tic
import Core.calc_overlap;
[mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,names2,namesBar);
toc


for ii=1:length(baridx2)
    baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
end
% we want to create a nicer structure for all-to-all comparison and
% contains easily accesible data.
tic
import Core.mp_res_to_struct;
[overlapStruct2] = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
toc

% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 2
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,xId,yId);
import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, overlapStruct2,xId,yId);


fig=figure
import CBT.Hca.UI.Helper.plot_best_pos;
plot_best_pos(fig, comparisonStruct([xId yId]), [], [], [],theoryStruct{1}.length,1/500*10^6);
% fig=figure
% import CBT.Hca.UI.Helper.plot_best_pos;
% plot_best_pos(fig, compStr([xId yId]), [], [], [],theoryStruct{1}.length,1/500*10^6);



maxpcc = PCC_OVERLAP(~isnan(PCC_OVERLAP));
import Zeromodel.beta_ev_params;
[parameters] = beta_ev_params(maxpcc, minOverlap/3);


pval = [];
for idxPair = 1:1000
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    import Zeromodel.beta_ev_cdf;
    pval(idxPair) = 1-beta_ev_cdf(overlapStruct(xId,yId).score, parameters(1)/MIN_OVERLAP_PIXELS * MIN_OVERLAP_PIXELS, 1, parameters(2),1);
    pval(idxPair)
end
% import Core.compute_dense_zvals
%   zVal(i,j) = compute_dense_zvals(overlapStruct(xId,yId).score,...
%  overlapStruct(xId,yId).lenA, overlapStruct(xId,yId).lenB, overlapStruct(xId,yId).h, 1, 'minimum')
% tic
% zVal = [];
% for i=1:size(overlapStruct,1)
%     i
%     for j=1:size(overlapStruct,2)
%           zVal(i,j) = compute_dense_zvals(overlapStruct(i,j).score,...
%  overlapStruct(i,j).lenA, overlapStruct(i,j).lenB, overlapStruct(i,j).h, 1, 'minimum');
%     end
% end
% toc

% import Core.plot_match_pcc;
% [f] = plot_match_pcc(barStruct, overlapStruct2,xId,yId);

%%
nT = 100;
tp = zeros(1,nT);
fp = zeros(1,nT);
tp2 = zeros(1,nT);
fp2 = zeros(1,nT);
fpList =[];fpList2=[];
tpList =[];tpList2=[];

for ii=1:nT
    idxPair = ii;
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    px1 = comparisonStruct{xId}.pos(1):comparisonStruct{xId}.pos(1)+comparisonStruct{xId}.lengthMatch-1;
    px2 = comparisonStruct{yId}.pos(1):comparisonStruct{yId}.pos(1)+comparisonStruct{yId}.lengthMatch-1;
    overLap = intersect(mod(px1,theoryStruct{1}.length),mod(px2,theoryStruct{1}.length));
    if length(overLap)>=100
        tp(ii)=1;
        tpList = [tpList xId yId];

    else
        fp(ii)=1;
        fpList = [fpList xId yId];
    end
% 
%         px1 = compStr{xId}.pB(1)-compStr{xId}.pA(1)+1:compStr{xId}.pB(1)-compStr{xId}.pA(1)+1+compStr{xId}.lengthMatch-1;
%         px2 = compStr{yId}.pB(1)-compStr{yId}.pA(1)+1:compStr{yId}.pB(1)-compStr{yId}.pA(1)+1+compStr{yId}.lengthMatch-1;
%     overLap = intersect(mod(px1,theoryStruct{1}.length),mod(px2,theoryStruct{1}.length));
%     if length(overLap)>=200
%         tp2(ii)=1;
%         tpList2 = [tpList2 xId yId];
% 
%     else
%         fp2(ii)=1;
%         fpList2 = [fpList2 xId yId];
% 
%     end
end

% f=figure,plot(cumsum(fp2))
% hold on
f=figure
plot(cumsum(fp))
% Unrecognized function or variable 'tp1'.

legend({'MP','PCC'})
saveas(f,'figs/fig13.png')

%% PLOT THRY COMP
idx = yId;
% imresize(barcodeGen{idx}.rawBarcode,'Scale' ,[0 comparisonStruct{idx}.bestBarStretch])
curBar = imresize(barcodeGen{idx}.rawBarcode(barcodeGen{idx}.rawBitmask),'Scale' ,[1 comparisonStruct{idx}.bestBarStretch]) 

if comparisonStruct{idx}.or(1)==2
    curBar = fliplr(curBar);
end

thr = importdata(theoryStruct{1}.filename);
% thrLambda = importdata(theoryStruct2{1}.filename);

% rawBg = barcodeGen{idx}.rawBg; %mode(cellfun(@(x) x.rawBg,barcodeGenLambda));
% meanlambda = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGenLambda));
% meanbar = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGen));

figure,plot(comparisonStruct{idx}.pos:comparisonStruct{idx}.pos+length(curBar)-1,zscore(curBar))
hold on
plot(comparisonStruct{idx}.pos:comparisonStruct{idx}.pos+length(curBar)-1,zscore(thr(comparisonStruct{idx}.pos:comparisonStruct{idx}.pos+length(curBar)-1)))

