% furst test for synthetic data
testSet={};
synth = 1;
if synth
    minLen = 2000;
    numF = minLen;
end

sF =1;
minOverlap = 300;%*0.95;

% sF = 1;
% xx = 300:10:1000;
xx = minLen;
psf =  300/110;

numEx = 10000;
import Zeromodel.prep_data;
[barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth,psf,numEx);

%  minO
xcorrs =[];
numElts = [];

for i=1:length(barcodeGen1)
    import SignalRegistration.masked_pcc_corr;
    [ xcorrs{i}, numElts{i} ] = masked_pcc_corr( barcodeGen1{i}.rawBarcode, barcodeGen2{i}.rawBarcode,...
       barcodeGen1{i}.rawBitmask,  barcodeGen2{i}.rawBitmask,minOverlap ); %todo: include division to k to reduce mem
%      xcorrs{i} =  xcorrs{i}.*numElts{i}./(numElts{j}-1); % change to N-1 facotr

end

import Zeromodel.beta_ev_fit;
import Zeromodel.beta_ev_nu;

% res = cell(1,length(barcodeGen1));
pars = zeros(size(xcorrs{1},2),1);
for j=1:size(xcorrs{1},2)
    j
%          xcorrs{i} =  xcorrs{j}.*numElts{j}./(numElts{j}-1);

    scores =  cellfun(@(x) x(1,j),xcorrs);
%     if length(scores) > 50

        [pars(j,1)] = beta_ev_nu(scores, minLen/3);
%         [pars(j,1),pars(j,2)] = beta_ev_fit(scores,[4 1],[inf inf],[4 1],[0 1]);
end

figure,histogram(pars(:,1),20)
[1/(sqrt(2*pi)*psf) mean(pars(:,1)/minLen)]

% todo: test based on different psf
%%

testSet ={ '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
synth = 0;
minLen = 400;
numF = 10;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

minOverlap = 300;
import Zeromodel.prep_data;
[barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

sF = 1;
tic
import Core.calc_overlap_pcc_sort;
[overlapStruct] = calc_overlap_pcc_sort([barcodeGen1 barcodeGen2], sF,minOverlap);
%%
overlapStruct1 = overlapStruct(end/2+1:end,1:end/2);
% overlapStruct1 = overlapStruct(1:length(barcodeGen1),1:length(barcodeGen1));

% overlapStruct1 = overlapStruct;

PCC_OVERLAP = reshape([overlapStruct1.score], size(overlapStruct1,1),size(overlapStruct1,2));
overlaplen = reshape([overlapStruct1.overlaplen], size(overlapStruct1,1),size(overlapStruct1,2));
bestBarStretch = reshape([overlapStruct1.bestBarStretch], size(overlapStruct1,1),size(overlapStruct1,2));

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen1,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen1,'un',false)]',{'rawBarcode','rawBitmask'},2);

 barStruct2 = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen2,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen2,'un',false)]',{'rawBarcode','rawBitmask'},2);


idxPair = 5;
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
% normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(filtM(:),'desc','MissingPlacement','last');

[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

% length(sortedVals(sortedVals>sortedVals2(2)))
% number of false positives based on two dataset, can then compare with
% pvals


import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct1,xId, yId,barStruct);
%%

% indScores = reshape({overlapStruct1.indScores}, size(overlapStruct1,1),size(overlapStruct1,2));
% indOverlaps = reshape({overlapStruct1.indOverlaps}, size(overlapStruct1,1),size(overlapStruct1,2));
% 
allscores = {overlapStruct1.score};
% lenBasedScore = arrayfun(@(x) cellfun(@(


tic
lenBasedScore = arrayfun(@(x) cellfun(@(y) y(x),allscores) , xx,'un',false);
toc
indScores2 = reshape([lenBasedScore{1}], size(overlapStruct1,1),size(overlapStruct1,2));

% indScores = reshape({overlapStruct1.indScores(xx)}, size(overlapStruct1,1),size(overlapStruct1,2));

%
import Zeromodel.beta_ev_params;

% % fit the parameters:
% maxpcc = PCC_OVERLAP(~isnan(PCC_OVERLAP));
% [parameters] = beta_ev_params(maxpcc, minOverlap);
% 
% parameters(1)/minOverlap
% parameters(2)
% 2*(RAND_LENGTH_MIN+RAND_LENGTH_MIN-2*minOverlap)

%
% overlap = overlaplen(~isnan(PCC_OVERLAP));

% pd = fitdist(overlap(:)-minOverlap+1,'exponential');


parameters1 =[];
parameters2 =[];

import Zeromodel.beta_ev_params;
import Zeromodel.beta_ev_fit;
import Zeromodel.beta_ev_fit_test;

% [a_fit, n_fit] = beta_ev_fit(scores)

res = cell(1,length(lenBasedScore));
for j=1:length(lenBasedScore)
    j
    scores = lenBasedScore{j}(~isnan(lenBasedScore{j}));
%     if length(scores) > 50

%         [parameters] = beta_ev_params(scores, minOverlap/3);
        [parameters(1), parameters(2)] = beta_ev_fit(scores);

        parameters1 = [parameters1 parameters(1)];
        parameters2 = [parameters2 parameters(2)];

%         m = bootstrp(10,@beta_ev_params,scores,minOverlap/3);
%         [m] = bootstrp(100,@beta_ev_fit_test,scores);
% 
%         res{j} = m;
end
%        

%% check dist for non matching real
overlapStruct1 = overlapStruct(end/2+1:end,1:end/2);

rgg = 300:100:1200;
ppvec =[];
for numPts = rgg;
    allpairs = overlapStruct1(:);
    sc = [];
    for i=1:10000
        newpts =allpairs(i).xcorrs(allpairs(i).numElts{1}==numPts);
        if length(newpts(:))==4
            sc =[ sc; newpts(:) ];
        end
    end
    import Zeromodel.beta_ev_nu;
    [pars] = beta_ev_nu(sc, numPts/3);
    
    ppvec = [ppvec pars/numPts];
    pars/numPts
end

% figure,plot(rgg,ppvec.*(300:10:800))
figure,plot(rgg,ppvec)

predPSF =1/(mean(ppvec)*sqrt(2*pi))

%% compare to thry

fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)

psffac = 0.5;
nmbp = barcodeGen1{1}.nmbp*psffac;
nmpx = 110*psffac;
import Thry.gen_theoretical;
[theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);

[comparisonStruct] = compare_to_t(barcodeGen1,theoryStruct,sF,sets)
mean(cellfun(@(x) x.maxcoef(1),comparisonStruct))
