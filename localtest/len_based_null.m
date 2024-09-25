testSet={};
synth = 1;
if synth
    minLen = 900;
    numF = 11000;
end

sF =1;
minOverlap = minLen;%*0.95;

% sF = 1;
% xx = 300:10:1000;
xx = minLen;
psf = 300/110;

import Zeromodel.prep_data;
[barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth,psf);

 



tic
import Core.calc_overlap_pcc_null;
[overlapStruct] = calc_overlap_pcc_null([barcodeGen1 barcodeGen2], sF,minOverlap);

overlapStruct1 = overlapStruct(end/2+1:end,1:end/2);
% overlapStruct1 = overlapStruct;

indScores = reshape({overlapStruct1.indScores}, size(overlapStruct1,1),size(overlapStruct1,2));
indOverlaps = reshape({overlapStruct1.indOverlaps}, size(overlapStruct1,1),size(overlapStruct1,2));
% 
allscores = {overlapStruct1.indScores};
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
%             [parameters] = beta_ev_params(maxpcc, minOverlap/3);

%     parameters = [xx/(sqrt(2*pi)*psf) max(4,indOverlaps{1,2}(xx)/(sqrt(2*pi)*psf))]
    xx1 =-0.3:0.01:1;
    import Zeromodel.beta_ev_pdf
    values = beta_ev_pdf(xx1, parameters(1), 1, parameters(2));
    figure,plot(xx1,values)

hold on
histogram(scores,'Normalization','pdf')
% legend({strcat(['EV dist fit, \nu=' num2str(parameters(1)) ' \lambda=' num2str(parameters(2))])},'location','southoutside')
% legend({strcat(['EV dist fit, overlap=' num2str(minOverlap) ', R_1=' num2str(parameters(1)./xx) ' R_2=' num2str(parameters(2)/indOverlaps{1,2}(xx))])},'location','southoutside')
legend({strcat(['EV dist fit, overlap=' num2str(minOverlap) ', R_1=' num2str(parameters(1)./xx) ' R_2=' num2str(parameters(2))])},'location','southoutside')

[1/(sqrt(2*pi)*psf) max(4,indOverlaps{1,2}(xx)/(sqrt(2*pi)*psf))]
% set(gca, 'YScale', 'log')
% end
%%
figure,
tiledlayout(2,1) 
nexttile
plot((xx),parameters1./(xx))
xlabel('\nu/overlap')
nexttile
plot((xx),parameters2)
xlabel('\lambda')

%
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

%% Distribution
% two different cases: no correlation & correlated, which is described by
% number of attempts to fit.
% When nA = 4, we estimate

import Zeromodel.beta_ev_cdf;

lenOverlap = 300;

A1 = sqrt(1/2*pi);
sigma = 2.3;
nu = lenOverlap*A1/sigma; % where sigman = PSF (in pixels) 
nA = 40;
% and
pv =[];
for i=0:0.01:1
    pv = [pv (1-beta_ev_cdf(i, nu, 1, nA,0))];
end
figure,plot(0:0.01:1,pv)



% dependence on overlap length + num of attempts + psf

for overlap=300:10:1000;
    for numAttempts = 4:10:1000

    end
end

% interpolation for speedy conversion

par1 = cellfun(@(x) x(1),fullmat);
par2 = cellfun(@(x) x(2),fullmat);

par1 = 30:1000;
par2 = 4:10:400;

par1Fun = @(x,y,z) interp3(length1, length2,w, par1, x,y,z,'nearest');
par2Fun = @(x,y,z) interp3(length1, length2,w, par2, x,y,z,'nearest');

p = @(x,y,z) 1-(0.5+0.5*(1-betainc((x).^2,0.5,y/2-1,'upper'))).^z ;


%%  
% for individual barcodes, use PCC distribution analysis, comparison
% between two dist
f = @(x,n) (1-x.^2).^((n-4)/2)./beta(1/2,1/2*(n-2));

n = 0.1;
x=-1:0.01:1;


f = @(x,sigma) 1/2*(1+erf(x/(sqrt(2)*sigma)));

p = @(x,y) 1/2*(1+sign(x).*betainc((x).^2,0.5,y/2-1,'lower'));
% p = @(x,y) betainc((x).^2,0.5,y/2-1,'lower');

rang =  -0.99:0.0001:0.99;
nu=1000;
figure,plot(rang,p(rang,nu))
hold on
plot(rang,f(rang,1/sqrt(nu-2)))

% figure,plot
figure,plot(rang,p(rang,nu)-f(rang,1/sqrt(nu-2)))
red=arrayfun(@(x) sum(abs(p(rang,nu)-f(rang,1/sqrt(nu-x)))),0:0.01:3);
tic; v =p(rang,400);toc;
tic; v =f(rang,200);toc;

  tic;                scaledXcorrs = f(xcorrs,numElts{s} ); toc

  xcorrs(xcorrs<0) = nan;
  xcorrs(isnan(xcorrs)) = 0;
  numElts{s}(isnan( numElts{s})) = 100;
%  tic;                scaledXcorrs = p(xcorrs(:),numElts{s}(:) ); toc
%                 f = @(x,n) (1-x.^2).^((n-4)/2)./beta(1/2,1/2*(n-2));
%                 scaledXcorrs = xcorrs./log(numElts{s});

                f = @(x,n) 1-1/2*(1+erf(x/(sqrt(2)*n)));

            toc

% figure,plot(x,f(x,n))

p2 = @(x,y,z) (1/2*(1+sign(x).*betainc((x).^2,0.5,y/2-1,'lower'))).^z;
% tic; v =p(rang,400);toc;
rang =  -0.99:0.0001:0.99;
nu=1000;lambda=4;
figure,plot(rang,p2(rang,nu,lambda))
