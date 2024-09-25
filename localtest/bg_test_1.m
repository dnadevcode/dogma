
% testSet ={'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
testSet = {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/'};
% testSet ={ '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
%     '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
synth = 0;
minLen = 200;
numF = 10;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

minOverlap = 150;
import Zeromodel.prep_data;
[barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

sF = 0.8:0.01:1.2;
tic
import Core.calc_overlap_pcc_sort_m;
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen1 barcodeGen2], sF,minOverlap);
% synth = 0;
% minLen = 400;
% numF = 10;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% 
% minOverlap = 300;
% import Zeromodel.prep_data;
% [barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)
% 
% sF = 1;
% tic
% import Core.calc_overlap_pcc_sort;
% [overlapStruct] = calc_overlap_pcc_sort([barcodeGen1 barcodeGen2], sF,minOverlap);
% 
% 
%%
overlapStruct1 = overlapStruct;

% % pvals
ccthresh = 0.7;
totL = length([overlapStruct1.sc]);
pscores = nan(1,totL);
nuF = 0.05; % dep on PSF.
import Zeromodel.beta_ev_cdf

for ii=1:totL
    ii
    if overlapStruct1(ii).sc > ccthresh
        pscores(ii) =-( 1-beta_ev_cdf(overlapStruct1(ii).sc, nuF*overlapStruct1(ii).overlaplen, 1, 2*max(overlapStruct1(ii).lenA,overlapStruct1(ii).lenB)));
    end
end

%%
% overlapStruct1 = overlapStruct(1:length(barcodeGen1),1:length(barcodeGen1));

% overlapStruct1 = overlapStruct;

% PCC_OVERLAP = reshape([overlapStruct1.sc], size(overlapStruct1,1),size(overlapStruct1,2));
PCC_OVERLAP = reshape([pscores], size(overlapStruct1,1),size(overlapStruct1,2));

overlaplen = reshape([overlapStruct1.overlaplen], size(overlapStruct1,1),size(overlapStruct1,2));
bestBarStretch = reshape([overlapStruct1.bestBarStretch], size(overlapStruct1,1),size(overlapStruct1,2));

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen1,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen1,'un',false)]',{'rawBarcode','rawBitmask'},2);

%  barStruct2 = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen2,'un',false);...
%     cellfun(@(x) x.rawBitmask,barcodeGen2,'un',false)]',{'rawBarcode','rawBitmask'},2);


idxPair = 1;
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
% normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(filtM(:),'desc','MissingPlacement','last');

[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

% length(sortedVals(sortedVals>sortedVals2(2)))
% number of false positives based on two dataset, can then compare with
% pvals


import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct1,xId, yId,barStruct,sortedVals(idxPair));
sortedVals(idxPair)
%%


fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)

psffac = 1;
nmbp = barcodeGen1{1}.nmbp*psffac;
nmpx = 110*psffac;
import Thry.gen_theoretical;
[theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);

sets.comparisonMethod = 'mass_pcc';

[comparisonStruct] = compare_to_t(barcodeGen1,theoryStruct,sF,sets);


idxPair = 4;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct1,xId, yId,barStruct);

bpPx = nmpx/nmbp;
 import CBT.Hca.UI.Helper.plot_best_pos;
    f = plot_best_pos([], comparisonStruct([xId yId]), [], [], [],theoryStruct{1}.length,1/1*10^6);
%  s

[comparisonStructMP] = compare_to_t_mp(barcodeGen1([xId,yId]),theoryStruct,sF,150);

 import CBT.Hca.UI.Helper.plot_best_pos;
    f = plot_best_pos([],comparisonStructMP, [], [], [],theoryStruct{1}.length,1/1*10^6);
%  
lenB = comparisonStructMP{1}.lengthMatch;

ix= 2 ;
N=min(150,length( comparisonStructMP{ix}.allposB));
pos = comparisonStructMP{ix}.allposB(1:N)-comparisonStructMP{ix}.allposA(1:N);
figure
plot(pos,'x')

%% the sample plot but barcode vs barcode
MIN_OVERLAP_PIXELS = 150;
mkdir(strcat('output',timestamp))
[overlapStruct2] = calc_overlap_mp(barcodeGen1([xId,yId]),overlapStruct1(xId,yId).bestBarStretch, MIN_OVERLAP_PIXELS,timestamp);

% [mpI1individual,mp1individual,maxMP,stridxindividual,compStr] = compare_mp_all(theoryStruct,barcodeGen,minLen,ix,timestamp,sF,MIN_OVERLAP_PIXELS,numWorkers)


ix= 2 ;
N=min(150,length( overlapStruct2(1,2).allposB));
pos =  overlapStruct2(1,2).allposB(1:N)-overlapStruct2(1,2).allposA(1:N);
figure
plot(pos,'x')
% overlapStruct1(xId,yId).sc
% comparisonStructMP{1}.allposA
%%
MIN_OVERLAP_PIXELS = 150;

N = 150;
NN = 2000;
curStd = zeros(1,NN);
for idxPair = 1:NN;
    idxPair
    
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    
    [overlapStruct2] = calc_overlap_mp(barcodeGen1([xId,yId]),overlapStruct1(xId,yId).bestBarStretch, MIN_OVERLAP_PIXELS,timestamp);
    curPts =  overlapStruct2(1,2).allposB(1:min(end,N))-overlapStruct2(1,2).allposA(1:min(end,N));
    curStd(idxPair)  = std( overlapStruct2(1,2).allposB(1:min(end,N))-overlapStruct2(1,2).allposA(1:min(end,N)));

end

thrStd = 50;
f=figure
plot(cumsum(curStd>thrStd))
%%
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
f=figure
% plot(cumsum(fp))
plot(fp)
%%
badpts = find(fp);

i=1;
    idxPair = badpts(i)
    
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct1,xId, yId,barStruct);
bpPx = nmpx/nmbp;
 import CBT.Hca.UI.Helper.plot_best_pos;
    f = plot_best_pos([], comparisonStruct([xId yId]), [], [], [],theoryStruct{1}.length,1/1*10^6);
%  s
    [overlapStruct2] = calc_overlap_mp(barcodeGen1([xId,yId]),overlapStruct1(xId,yId).bestBarStretch, MIN_OVERLAP_PIXELS,timestamp);
    curPts =  overlapStruct2(1,2).allposB(1:min(end,N))-overlapStruct2(1,2).allposA(1:min(end,N));

    %%