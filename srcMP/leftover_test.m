% this will contain the leftover_score test

%{
Run leftover test. 
1) Calculate all overlaps.
2) Extract best match
3) Bootstrap partial overlap score for best matches (and different lengths)
4) Fit EVD for each partial overlap score histogram
5) Convert actual partial scores to p-values
6) Do the same for local scores.


idea: for partial scores, we also bootstrap, to the same length as the
local overlap. This works as long as we have enough pixels in leftover
vector (more than 10/50?)
%}
s = RandStream('mlfg6331_64'); 

idxRun = 5;
oS= resRun{idxRun}.oS;
barcodeGen = bG{idxRun};

[localPCC, sortedIdsAll, pscores,fullPCC,overlaplen,partialPCC,partialLength] = sorted_scores(oS);

N = 2000;
pairedBar = zeros(2,oS(1,2).h,N);
pairedPartialBar = cell(1,N);
pairedFullBar = cell(1,N);
for idxPair = 1:N;
    [curSink, curSource] = ind2sub(size(oS),sortedIdsAll(idxPair));
    % bootstrap
    k = curSink;
    iy = curSource;
    
    import Core.get_full_overlap_score;
    [fullscore,overlaplen,lenB, lenA, partialPCC,parLen,aM,bM,aP,bP] = get_full_overlap_score(oS(k,iy).pA,oS(k,iy).pB,...
             oS(k,iy).bestBarStretch, oS(k,iy).or,barcodeGen([k iy])',oS(k,iy).h);
    pairedBar(:,:,idxPair) = [aM;bM];
    pairedPartialBar{idxPair} = [aP;bP];
    pairedFullBar{idxPair} = [aM;bM];

end

% bootstrap
B = reshape(pairedBar, [2, oS(1,2).h*N]);


wBoostrapping = 200;
% now split into ~50px windows. skips the last window if not divisible
% numWindows = floor(length(B)/wBoostrapping);
    
%     R = reshape(1:floor(length(B)/numWindows)*numWindows,floor(length(B)/numWindows),[]);
     
    
NN = 100;
pccScore = zeros(1,NN);
for i=1:NN
    y = datasample(s,1:length(B),wBoostrapping,'Replace',true);

    b1 = B(1,y);
    b2 = B(2,y);
    pccScore(i) = zscore(b1(:)',1)*zscore(b2(:),1)/length(b1(:)); % todo: speed up by pre-calculating the residuals and only scale by mean

end
% figure,histogram(pccScore)

NN = 100;
wBoostrapping = 100;

import Nullmodel.leftover_bootstrapping;
minL = 50;
pccScoreB = leftover_bootstrapping(wBoostrapping,pairedPartialBar,partialLength, NN,minL);
% figure,histogram(cellfun(@(x) mean(x),pccScoreB))

scoresAll = cellfun(@(x) x(1),pccScoreB);
scores = scoresAll(300:end);
figure,histogram(scores,'Normalization','pdf')
hold on
pccScoreB = leftover_bootstrapping(500,pairedPartialBar,partialLength, NN,minL);
scoresAll = cellfun(@(x) x(1),pccScoreB);
scores = scoresAll(300:end);
histogram(scores,'Normalization','pdf')

scores = scores(~isnan(scores));


import Nullmodel.full_dist_mle;
f = @(x0) full_dist_mle(scores(:),0,x0);

 nFit = [fsolve(f,wBoostrapping/4,optimoptions('fsolve','Display','off'))]

 import Nullmodel.full_dist;

 ccv = -0.99:0.01:0.99;
[fvals] = arrayfun(@(x) full_dist(x, 0, nFit),ccv)
% figure,plot(ccv,fvals)
figure,histogram(scoresAll,'Normalization','pdf')
hold on
plot(ccv,fvals,'red','LineWidth',2)
legend({'Left-over dist','fit'})

 %
 %            parameters = [fsolve(f,x0,optimoptions('fsolve','Display','off'))];

% do similar function for global

%         b2 = bar2(R(:,y));
%         pccScore(i) = zscore(b1(:)',1)*zscore(b2(:),1)/length(b1(:)); % todo: speed up by pre-calculating the residuals and only scale by mean
    
%% test to  see dependence on n_effectve of variance
import Nullmodel.full_dist_mle;
import Nullmodel.leftover_bootstrapping;
import CBT.Hca.Core.Pvalue.compute_evd_params;
NN = 100;
wBoostrapping = [10:50:500];
stdW = zeros(1,length(wBoostrapping));
nFit = zeros(1,length(wBoostrapping));
parameters = zeros(length(wBoostrapping),2);
minL = 50;
nSt = 500;
h=oS(1,2).h;

for  i= 1:length(wBoostrapping);
    i
    pccScoreB = leftover_bootstrapping(wBoostrapping(i),pairedPartialBar(nSt:end),partialLength(nSt:end), NN,minL);
    
    scoresAll = cellfun(@(x) x(1),pccScoreB);
    scoresAll = scoresAll(~isnan(scoresAll));
    % scores = scoresAll(300:end);
    % figure;histogram(scoresAll,'Normalization','pdf')
    stdW(i) = nanstd(scoresAll);
    
    % for leftovers
    f = @(x0) full_dist_mle(scoresAll(:),0,x0);
    nFit(i) = [fsolve(f,wBoostrapping(i)/4,optimoptions('fsolve','Display','off'))];
    % for full
    pccScoreFull = leftover_bootstrapping(wBoostrapping(i),pairedFullBar(nSt:end),h*ones(length(pairedFullBar(nSt:end)),1), NN,minL);

    scoresAllFull = cellfun(@(x) x(1),pccScoreFull);
    scoresAllFull = scoresAllFull(~isnan(scoresAllFull));
 
     [ parameters(i,:) ] = compute_evd_params( scoresAllFull,wBoostrapping(i)/4 );
end

figure,plot(wBoostrapping,nFit)
figure,plot(wBoostrapping,stdW)
hold on
plot(wBoostrapping,1./sqrt(nFit))
plot(wBoostrapping,1./sqrt(parameters(:,1)))
xlabel('Bootstrapping window')
legend({'Bootstrapped st.dev.','1./sqrt(nFit)','1/sqrt(nu)'})
% figure,plot(wBoostrapping,stdW)

scores = scores(~isnan(scores));


%  nFit = [fsolve(f,wBoostrapping/4,optimoptions('fsolve','Display','off'))]

pdfF = @(cc,evdPar) evdPar(2)*(1/2*(1+betainc(cc.^2,1/2, evdPar(1)/2-1))).^(evdPar(2)-1).*(1-cc.^2).^((evdPar(1)-4)/2)./beta(1/2,  evdPar(1)/2-1);





[histAll,vals] = histcounts(scoresAllFull,'Normalization','count');
bincenters = (vals(2:end)+vals(1:end-1))/2;
cc = bincenters;
pdfPlot = arrayfun(@(x) pdfF(x, [parameters(i,1)  parameters(i,2)]),cc);

% implement chi^2 test for the goodness-of-fit
import Zeromodel.chi2_test;
[chi2Score] = chi2_test(histAll, pdfPlot,[1 length(histAll)],5); % make a histogram of rand scores?

figure,histogram(scoresAllFull,'Normalization','pdf')
hold on
plot(cc,pdfPlot)

