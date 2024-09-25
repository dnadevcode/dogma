% 09/06/22 True positive test

testSet = {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220221_Sample-dEC-st10_464.93bpPERpx_0.236nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220224_Sample-dEC-st11_425.53bpPERpx_0.258nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220223_Sample-dEC-st11_519.34bpPERpx_0.211nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220225_Sample-dEC-st11_553bpPERpx_0.198nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220613_Sample DA32087-4-st2_570bpPERpx_0.193nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/Final_Data_Shrink_Sorted/20220610_Sample DA32087-4-st1_570bpPERpx_0.192nmPERbp/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/ch1_test/'};

idx = 7;
numF = 1;

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');


import Core.create_barcodegen;
[barcodeGen, kymoStructs, sets] = create_barcodegen(testSet{idx},1,numF,1,timestamp,1);
load('/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/barcodeGen993.mat')


% only keep barcodes longer than minimum length.
minLen = 200; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(barLens>minLen);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);


sets.comparisonMethod = 'mass_pcc';
sets.filter = 0;
sets.filterSettings.filter = 0;




%%
fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};
fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/sequence1.fasta',...
    '/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/sequence2.fasta',...
    '/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/sequence3.fasta'};

%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)

% nmbp = barcodeGen{1}.nmbp;
% nmbp = 0.225;
nmbp = 0.2508;
nmpx = 254;
psffac = 300/681.81;
nmbp = nmbp*psffac;
nmpx = nmpx*psffac;

import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmpx);

sum(cellfun(@(x) x.length,theoryStruct))
%% PCC /can compare sliding instead of hca
sF = 0.9:0.01:1.1;
% sF = 1;
minOverlap = 100;
tic
import Core.calc_overlap_pcc_individual;
[oS,overlapStruct] = calc_overlap_pcc_individual(barcodeGen,barcodeGenT, sF,minOverlap);

bpPx = nmpx/nmbp;
import CBT.Hca.UI.Helper.plot_best_pos;
f = plot_best_pos([], overlapStruct(:,1)', [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%
PCC_OVERLAP = reshape([oS.sc], size(oS,1),size(oS,2));
overlaplen = reshape([oS.overlaplen], size(oS,1),size(oS,2));

filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');
%%
% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 10;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));


PCC_OVERLAP(xId,yId)

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

barStruct2 = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGenT,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGenT,'un',false)]',{'rawBarcode','rawBitmask'},2);
% have to save separately..

import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, oS,xId, yId,barStruct2);
hold on
xlim([oS(xId,yId).pA-100 oS(xId,yId).pA+oS(xId,yId).lenA+100])
%
%%
% nmbp = barcodeGen{1}.nmbp;
% nmbp = 0.225;
% nmpx = 170; % 208?
% psf = 300; %300 nm
% psffac = 1;
% nmbp = nmbp*psffac;
% nmpx = nmpx*psffac;
% 
% import Thry.gen_theoretical;
% [theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmpx);

sF = 0.9:0.01:1.1;
% sF = 1;

rezMax=[];bestBarStretch=[];bestLength=[]; % todo: write this as overlapstructure?
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

mean(cellfun(@(x) x.maxcoef(1),comparisonStruct))

sets.timeFramesNr = nan;
sets.displayResults=1
sets.userDefinedSeqCushion = 0;
sets.genConsensus = 0;
import CBT.Hca.UI.get_display_results;
[res] = get_display_results(barcodeGen,[], comparisonStruct, theoryStruct, sets);
%%
xId = 242;
yId= 3;
import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, oS,xId, yId,barStruct2);
hold on
xlim([oS(xId,yId).pA-100 oS(xId,yId).pA+oS(xId,yId).lenA+100])


%%
% %% MP
numWorkers = 30;
MIN_OVERLAP_PIXELS_BL = 300;
minLen_BL = 300;
% minLen = 500; %150kb? or less? depends on application // if GUI, user selects thi
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
[a,b] = find(barLens>minLen_BL);
% find(b==44)
import Core.compare_mp_all
ix = 1;
[mpI1individual,mp1individual,maxMP,stridxindividual,compStr] = ...
    compare_mp_all(theoryStruct,barcodeGen,minLen_BL,ix, timestamp,sF,MIN_OVERLAP_PIXELS_BL,numWorkers);

lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));
% lengthBorders = lengthBorders(ix);
% fig1 = figure;
bpPx = nmpx/nmbp;
% import CBT.Hca.UI.Helper.plot_best_pos;
% f= plot_best_pos([], compStr, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
% % 
% saveas(f,'figs/fig10.png')

% plot_best_pos([], comparisonStruct, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

%% Overlaps
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


numWorkers = 30;
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
[overlapStruct] = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
toc
% import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, overlapStruct,11,6);

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
idxPair = 1;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));


PCC_OVERLAP(xId,yId)

import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, overlapStruct,xId,yId);

% overlapStruct(xId,yId).pB-overlapStruct(xId,yId).pA+1
% 
% [res] = get_display_results(barcodeGen([xId yId]),[], comparisonStruct([xId yId]), theoryStruct, sets);
  import CBT.Hca.UI.Helper.plot_best_pos;
    f = plot_best_pos([], compStr([find(b==xId) find(b==yId)]), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%  saveas(f,'figs/fig12.png')
% 
    plot_best_pos([], comparisonStruct([xId yId]), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%     plot_best_pos([], comparisonStruct, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

% % % 
pthresh = 0.000001;
import Core.simple_bargroup;
[orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,yId,pthresh,1);
[orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,21,pthresh,1);

% % % [orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,7,pthresh,1);
% % 
% %     plot_best_pos([], comparisonStruct(bars), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%     plot_best_pos([], compStr(bars), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

%%
nT = 300;
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

        px1 = compStr{xId}.pB(1)-compStr{xId}.pA(1)+1:compStr{xId}.pB(1)-compStr{xId}.pA(1)+1+compStr{xId}.lengthMatch-1;
        px2 = compStr{yId}.pB(1)-compStr{yId}.pA(1)+1:compStr{yId}.pB(1)-compStr{yId}.pA(1)+1+compStr{yId}.lengthMatch-1;
    overLap = intersect(mod(px1,theoryStruct{1}.length),mod(px2,theoryStruct{1}.length));
    if length(overLap)>=200
        tp2(ii)=1;
        tpList2 = [tpList2 xId yId];

    else
        fp2(ii)=1;
        fpList2 = [fpList2 xId yId];

    end
end

f=figure,plot(cumsum(fp2))
hold on
plot(cumsum(fp))
% Unrecognized function or variable 'tp1'.

legend({'MP','PCC'})
saveas(f,'figs/fig13.png')

%% just plots
fig =figure;
import CBT.Hca.UI.Helper.plot_best_pos;
plot_best_pos(fig, comparisonStruct, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,comparisonStruct,theoryStruct,timestamp);


overlapPairs = unique(tpList);
[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen(overlapPairs),comparisonStruct(overlapPairs),theoryStruct,timestamp);

fig=figure
import CBT.Hca.UI.Helper.plot_best_pos;
plot_best_pos(fig, comparisonStruct(overlapPairs), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);


thresh = 0.001;
% 
% f=figure,plot(cumsum(pval(fpList2) < thresh))
% hold on
% plot(cumsum(pval(fpList) < thresh))
% % Unrecognized function or variable 'tp1'.

legend({'MP','PCC'})
saveas(f,'figs/fig14.png')
