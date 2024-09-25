
% Scipt for generating null mod distribution

testSet ={ '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
numF = 10;
testSet2 = {'C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast 280721\'};


timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

minLen = 500;
synth = 0;
if synth
    minLen = 1000;
    numF = 1000;
end

import Zeromodel.prep_data;
[barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

    
minOverlap = 300;
sF =1;

tic
import Core.calc_overlap_pcc_null;
[overlapStruct] = calc_overlap_pcc_null([barcodeGen1 barcodeGen2], sF,minOverlap);
toc
indScores = reshape({overlapStruct.indScores}, size(overlapStruct,1),size(overlapStruct,2));
indOverlaps = reshape({overlapStruct.indOverlaps}, size(overlapStruct,1),size(overlapStruct,2));
% 
allscores = {overlapStruct.indScores};
% lenBasedScore = arrayfun(@(x) cellfun(@(y) y(x),allscores) ,1:length(allscores{1}),'un',false);
tic
lenBasedScore = arrayfun(@(x) cellfun(@(y) y(x),allscores) , 300:310,'un',false);
toc
% indScores = cellfun(@(x) x.indScores,overlapStruct,'un',false);
% indOverlaps = cellfun(@(x) x.indOverlaps,overlapStruct,'un',false);
% 
% tic
% import Core.calc_overlap_pcc;
% [overlapStruct] = calc_overlap_pcc([barcodeGen1 barcodeGen2], sF,minOverlap);
% toc

overlapStruct1 = overlapStruct(1:length(barcodeGen1),1:length(barcodeGen1));
if synth
    overlapStruct1 = overlapStruct;
else
    overlapStruct1 = overlapStruct(length(barcodeGen1)+1:end,1:length(barcodeGen1));
end
%%
PCC_OVERLAP = reshape([overlapStruct1.score], size(overlapStruct1,1),size(overlapStruct1,2));
overlaplen = reshape([overlapStruct1.overlaplen], size(overlapStruct1,1),size(overlapStruct1,2));
bestBarStretch = reshape([overlapStruct1.bestBarStretch], size(overlapStruct1,1),size(overlapStruct1,2));

if synth
     barStruct = cell2struct([cellfun(@(x) x.rawBarcode,[barcodeGen1 barcodeGen2],'un',false);...
    cellfun(@(x) x.rawBitmask,[barcodeGen1 barcodeGen2],'un',false)]',{'rawBarcode','rawBitmask'},2);
    barStruct2 = barStruct;
else
 barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen1,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen1,'un',false)]',{'rawBarcode','rawBitmask'},2);

 barStruct2 = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen2,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen2,'un',false)]',{'rawBarcode','rawBitmask'},2);
end
%%
idxPair = 20000;
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
% normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length
 [sortedVals, sortedIds] = sort(filtM(:),'desc','MissingPlacement','last');

[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

% length(sortedVals(sortedVals>sortedVals2(2)))
% number of false positives based on two dataset, can then compare with
% pvals


import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct2, overlapStruct1,xId, yId,barStruct);
%%
saveas(f,'fig1n2.eps','epsc')

%%


% we can check the scores in the matrix representing only ECOLI1 vs ECOLI2
% idx = 4;

% 
% 
% import Zeromodel.gen_overlaps_for_barcodegen;
% 
% % % test: generate p-values for bargroups
% % if ispc
% %     SCAMP_LINE = '.\..\SCAMP\build\Release\SCAMP.exe'; % windows
% % else
% %     sets.SCAMP_LINE = '~/SCAMP/';
% %     SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
% % end
% 
% % delete(gcp('nocreate'))
% % numWorkers = 32;
% % parpool('local',numWorkers)
% %%
% NUM_RAND_FRAGMENTS = 500;
% PSF_WIDTH_PIXELS = 300/110;
% RAND_LENGTH_MIN = 1500;
% RAND_LENGTH_2 = 800;
% 
% % sF = 1;%
% 
% sF = 0.8:0.01:1.2;
% 
% import Nullmodel.gen_random;
% [~,barcodeGen] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,0);
% %
% 
% minOverlap = 300;
% tic
% import Core.calc_overlap_pcc;
% [overlapStruct] = calc_overlap_pcc(barcodeGen, sF,minOverlap);
% toc

%%
import Zeromodel.beta_ev_params;

% fit the parameters:
maxpcc = PCC_OVERLAP(~isnan(PCC_OVERLAP));
[parameters] = beta_ev_params(maxpcc, minOverlap);
% 
% parameters(1)/minOverlap
% parameters(2)
% 2*(RAND_LENGTH_MIN+RAND_LENGTH_MIN-2*minOverlap)

%
overlap = overlaplen(~isnan(PCC_OVERLAP));

pd = fitdist(overlap(:)-minOverlap+1,'exponential');


parameters1 =[];
parameters2 =[];

for j=minOverlap:800
    scores = maxpcc(overlap==j);
    if length(scores) > 50
        import Zeromodel.beta_ev_params;

        [parameters] = beta_ev_params(scores, minOverlap/3);
        parameters1 = [parameters1 parameters(1)];
        parameters2 = [parameters2 parameters(2)];
    else
        break
    end
end

x =minOverlap:minOverlap+length(parameters1)-1;
cnu = polyfit(x,parameters1,1);
y_est1 = polyval(cnu,x);
clambda = polyfit(x,parameters2,1);
y_est2 = polyval(clambda,x);
% f1 = figure,histogram(overlap(:),300)
% saveas(f1,'fig2.eps','epsc')

import Zeromodel.beta_ev_cdf;

N = size(overlaplen,1)*size(overlaplen,2);
pval1 = [];
pval2 = [];

for idxPair = 1:10000
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    p1 = polyval(cnu,overlaplen(xId,yId));
    p2 = polyval(clambda,overlaplen(xId,yId));
    scale = pdf(pd,overlaplen(xId,yId)-minOverlap+1)*N;
    pval2(idxPair) = (1-beta_ev_cdf(overlapStruct(xId,yId).score, p1, 1, p2,0));
    pval1(idxPair) = scale*pval2(idxPair);
end

%%
dist = 0.00001:0.00001:0.01;
vals= arrayfun(@(x) length(pval2(pval2<x)),dist);
figure,plot(dist,vals,'red')
hold on
plot(dist,dist*N,'black')
legend({'FP','Theoretical FP'})
% for thresh = dist
% thresh = 0.00001;
% thresh*N
% pd = fitdist(r,'exponential');
f=figure; tiledlayout(2,2,'TileSpacing','tight','Padding','tight')
nexttile
histogram(maxpcc,50);title('PCC histogram')
nexttile
h = histfit(overlap(:)-minOverlap+1,50,'exponential')
% xlabel('Length parameter')
title('Length histogram')

% saveas(f,'fig2e.eps','epsc')

%     f = figure,plot(x,parameters1)
% y
u=1;

nexttile
plot(x,parameters1./u)
hold on
plot(x,y_est1./u,'r--','LineWidth',2)
title('\nu parameter')

% saveas(f,'fig4e.eps','epsc')


nexttile
plot(x,parameters2./u)
hold on
plot(x,y_est2./u,'r--','LineWidth',2)
plot(x,2*(3000 - 2*x-1))
title('\lambda parameter')
xlabel('Overlap length')
saveas(f,'fig5s.eps','epsc')
% 2*

%%
overlapStructGood = overlapStruct(1:length(barcodeGen1),1:length(barcodeGen1));
PCC_OVERLAP_SAMPLE = reshape([overlapStructGood.score], size(overlapStructGood,1),size(overlapStructGood,2));
overlaplen_SAMPLE = reshape([overlapStructGood.overlaplen], size(overlapStructGood,1),size(overlapStructGood,2));

filtM = triu(PCC_OVERLAP_SAMPLE);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');


% import Zeromodel.beta_ev_cdf;
% pval(idxPair) = 1-beta_ev_cdf(overlapStruct(xId,yId).score, parameters(1)/MIN_OVERLAP_PIXELS * MIN_OVERLAP_PIXELS, 1, parameters(2),1);

import Zeromodel.beta_ev_cdf;

N = size(overlaplen_SAMPLE,1)*size(overlaplen_SAMPLE,2);
pval = [];
for idxPair = 1:1000
    [xId,yId] = ind2sub(size(PCC_OVERLAP_SAMPLE),sortedIds(idxPair));
     p1 = polyval(cnu,overlaplen_SAMPLE(xId,yId));
     p2 = polyval(clambda,overlaplen_SAMPLE(xId,yId));
     scale = pdf(pd,overlaplen_SAMPLE(xId,yId))*N;
    pval(idxPair) = scale*(1-beta_ev_cdf(overlapStructGood(xId,yId).score, p1, 1, p2,1));
    pval(idxPair)
end

f = figure
plot(pval(1:100))


idxPair = 2;
[xId,yId] = ind2sub(size(PCC_OVERLAP_SAMPLE),sortedIds(idxPair));

% length(sortedVals(sortedVals>sortedVals2(2)))
% number of false positives based on two dataset, can then compare with
% pvals


import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct1,xId, yId,barStruct,pval(idxPair));
saveas(f,'figpvalb.eps','epsc')

%%
