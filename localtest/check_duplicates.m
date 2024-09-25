
% testSet ={'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
testSet = {'/proj/snic2022-5-384/users/x_albdv/data/CHR/Mapping_Old E.coli_ EF365_after shrink finder/'};
% testSet ={ '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
%     '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
synth = 0;
minLen = 100;
numF = 1;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

minOverlap = 100;
import Zeromodel.prep_data;
[barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

sF = 1;
tic
import Core.calc_overlap_pcc_sort_m;
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen1 barcodeGen2], sF,minOverlap);
% synth = 0;
%

%%
PCC_OVERLAP = reshape([overlapStruct.sc], size(overlapStruct,1),size(overlapStruct,2));

% overlaplen = reshape([overlapStruct1.overlaplen], size(overlapStruct1,1),size(overlapStruct1,2));
% bestBarStretch = reshape([overlapStruct1.bestBarStretch], size(overlapStruct1,1),size(overlapStruct1,2));

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen1,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen1,'un',false)]',{'rawBarcode','rawBitmask'},2);

%  barStruct2 = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen2,'un',false);...
%     cellfun(@(x) x.rawBitmask,barcodeGen2,'un',false)]',{'rawBarcode','rawBitmask'},2);


idxPair = 3;
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
% normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(filtM(:),'desc','MissingPlacement','last');

[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

% length(sortedVals(sortedVals>sortedVals2(2)))
% number of false positives based on two dataset, can then compare with
% pvals


import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,xId, yId,barStruct,sortedVals(idxPair));