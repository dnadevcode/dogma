%% pval_test_synth_local

% 1) Generate synth data
rng("default")

import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 100;
% minOverlapPx = 300;
PSF_WIDTH_PIXELS = 300/110;
randLenPx = 5000;
randLenPx2 = 5000;

[~ , barSynth] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,randLenPx2,0);
[~ , barSynth2] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,randLenPx,0);

% 2) Gen scores
minOverlapLen = 300;
randLenMin = 1600;
randLen2 = 2000;
sF = 1; %0.95:0.01:1.05;%

barLens = cellfun(@(x) sum(x.rawBitmask),barSynth);
barcodeGenGood = barSynth(barLens>randLenMin);
barLens2 = cellfun(@(x) sum(x.rawBitmask),barSynth2);
barcodeGenGood2 = barSynth2(barLens2>randLenMin);

%  calculates all MP and MPI
NN = 30; % NN changes the length randLen, 400:10:400+(NN-1)*10, so this changes the total number of attempts, but keeps same overlap length,
TT = 10; % changes in w
% so only n_fit should change
gap = 100;
out='output';
tic
import Nullmodel.gen_scores_real;
[maxPCC,totLen] = gen_scores_real(barcodeGenGood,barcodeGenGood2,randLenMin,randLen2,minOverlapLen,NN,sF,out,gap, TT,'SCAMP');
toc % todo: convert from local to full overlap!

totLenLeft =  (randLenMin:10:randLenMin+(NN-1)*10)-minOverlapLen+1;
overlapL = minOverlapLen+gap*([1:TT]-1);


save('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S2_data.mat','maxPCC','totLen','randLenMin','randLen2','minOverlapLen','NN','sF','out','gap', 'TT');
%% 

import Zeromodel.beta_ev_params;

casetest=1;
a_fit2 = zeros(NN,length(maxPCC{1}));
n_fit2 = zeros(NN,length(maxPCC{1}));
for k=1:NN
    for ii=1:length(maxPCC{k})


        if casetest==1
            [parameters] = beta_ev_params(maxPCC{k}{ii},(minOverlapLen+(ii-1)*gap ) /2);%(randLenMin+(k-1)*gap)/2);
%         else
%             N_FIXED = totLen(k)*(randLenMin+(k-1)*gap-minOverlapLen+1);
%             [parameters] = beta_ev_params(maxPCC{k}{ii}, 200, N_FIXED);
        end
        a_fit2(k,ii) = parameters(1);
        n_fit2(k,ii) = parameters(2);
        %     [a_fit(k), n_fit(k)] = beta_ev_fit(maxPCC{k}, [2 1], [inf inf], [4 (RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS)], [true false]);
    end
end

a_fit2 = a_fit2';
n_fit2 = n_fit2';

%% Similar plot for local pcc
barL2 = randLen2;
barL1 = (randLenMin+([1:NN]-1)*gap);

% figure,plot(mean(n_fit2')./(shortL-minOverlapLen))

par1 = zeros(1,TT);
par2 = zeros(1,TT);
par1CI = zeros(1,TT);
% choose 1
for j = 1:TT
%     aPar = a_fit2(j,:)/minOverlapLen;
    kk=1:NN;
    nPar = n_fit2(j,:);
    parEq = (2*(RAND_LENGTH_2+gap*(kk-1)));
    
    % fnpar = polyfit(parEq,nPar,1) % fit a line.    
    mdl = fitlm(parEq,nPar,'intercept',false);
    fnpar(1) = mdl.Coefficients.Estimate(1);
    fnpar(2) = 0;

    par2(j) = [mdl.Coefficients.Estimate(1)];
    par1(j) = [mean(a_fit2(j,:))]; % par1 mean
    par1CI(j) =  1.98*[std(a_fit2(j,:))/sqrt(length(a_fit2(j,:)))];
end

% par3 =  zeros(1,NN);
% for j=1:NN
%     nPar = n_fit2(:,j);
%     par3(j) = mean(nPar);
% end
%%
ix = 8;
par3 = n_fit2(ix,:); % can only do par3 for fixed w or fixed len2, since it depends on both. Write this in the text!

% minL = randLenMin+gap*([1:TT]-1);

betanu = 0.14;
betaN = 0.001;
f=figure('Position',[1 1 800 800])
tiledlayout(3,2,'TileSpacing','tight')

nexttile
hold on
plot(overlapL,par1./overlapL)

% errorbar(overlapL,par1./overlapL,par1CI./overlapL)
xlabel('$w$','Interpreter','latex')
ylabel('$\nu_{eff}/w$','Interpreter','latex')
title('(A) $\beta_{\nu}$ estimation','Interpreter','latex')
plot(overlapL,betanu*ones(1,length(overlapL)),'--')

nexttile
hold on
% plot(barL1,par3./(2.*(max(barL1-overlapL(ix),barL2-overlapL(ix)))))
plot(barL1,par3./((2.*((barL1-overlapL(ix)).*((barL2-overlapL(ix)))))))
% plot(barL1,par3./((2.*((barL1).*((barL2))))))

xlabel('$L_A$','Interpreter','latex')
ylabel('$N_{eff}/(2(\cdot L_B-w)(L_A-w)$','Interpreter','latex') %(2(\cdot L_B-w)(L_A-w)
title('(B) $\beta_{N}$ estimation','Interpreter','latex')
% 
plot(barL1,betaN*ones(1,length(barL1)),'--')

nexttile
% shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% longL = RAND_LENGTH_2+(0:NN-1)*gap;
imagesc(overlapL,barL1,a_fit2);
title('(C) $\nu_{eff}$ heatmap','Interpreter','latex')
colorbar;colormap(gray) 
ylabel('$L_A$','Interpreter','latex')
xlabel('$w$','Interpreter','latex')
nexttile
% shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% longL = RAND_LENGTH_2+(0:NN-1)*gap;
imagesc(overlapL,barL1,n_fit2);
colorbar;colormap(gray) 
title('(D) $N_{eff}$ heatmap','Interpreter','latex')
ylabel('$L_A$','Interpreter','latex')
xlabel('$w$','Interpreter','latex')

nexttile
% figure
iy = 30;
j = 1;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
% f=figure,
plot(xx,p)
hold on
histogram(maxPCC{iy}{j},'Normalization','pdf')
xlabel('Max local PCC','Interpreter','latex')
title('E) EVD fit histogram','Interpreter','latex')
xlim([0.5 0.8])
lgd = legend({'Functional EVD','PCC histogram'},'Interpreter','latex','Location','southoutside')


print('FIGS/FigS2.eps','-depsc','-r300');

%%
runextra = 0;
if runextra
%%%%%%%% 

% scNeff = 0.03;
scNeff = mean(a_fit2(:)./minOverlapLen);
nPar = mean(n_fit2,2);
parEq = 2*totLenLeft;
% parEq = 2*totLenLeft*(randLen2-minOverlapLen+1);
% parEq = 2*(totLenLeft+randLen2-minOverlapLen+1);%+2*scNeff*(randLen2-minOverlapLen+1);

% fnpar = polyfit(parEq,nPar,1) % fit a line.
mdl = fitlm(parEq,nPar,'intercept',false);
fnpar(1) = mdl.Coefficients.Estimate(1);
fnpar(2) = 0;


f=figure,
tiledlayout(2,2)
nexttile
plot(parEq,nPar)
hold on
title('A) EVD fit real data $N_{eff}$','Interpreter','latex')

kk=1:NN;
plot(parEq, fnpar(1)*parEq+fnpar(2))
xlabel('$N$','Interpreter','latex')
ylabel('$N_{eff}$','Interpreter','latex')
lgd = legend({'$N_{eff}$',['Line f(x) = ',num2str(fnpar(2)),' + ', num2str(fnpar(1)), 'x']},'Interpreter','latex')
lgd.Location = 'southoutside';
nexttile
aPar = mean(a_fit2,2)/minOverlapLen;
c1 = polyfit(1:length(aPar),aPar,1) % fit a line.
coef = aPar/minOverlapLen
plot(parEq,aPar)
hold on
plot(parEq,c1(1)*(1:length(aPar)) +c1(2))
lgd2 = legend({'$C_{\nu_{eff}}$',['Line f(x) = ',num2str(c1(2)),' + ', num2str(c1(1)), 'x']},'Interpreter','latex')
lgd2.Location = 'southoutside';
title('B) EVD fit real data $\nu_{eff}$','Interpreter','latex')


nexttile

iy=1;
ix=2;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(iy,ix), 1, n_fit2(iy,ix));
% f=figure,
plot(xx,p)
hold on
histogram(maxPCC{iy}{ix},'Normalization','pdf')
title('EVD fit histogram','Interpreter','latex')
% saveas(f,fullfile(foldsave,'figT2.png'))

% set(gca,'YScale','log')

%%


%%%%%%%%%%%%%%%%%%%

% Previous
%{
addpath(genpath(pwd));
% test: generate p-values for barcodes
SCAMP_LINE = '~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
delete(gcp('nocreate'))
numWorkers = 24;

parpool('local',numWorkers)
% This will depend on these parameters: overlap length, number of attempts, re-scaling factors, psf.
% possibly re-scaling factors and psf have minor relevance, then it would
% be mainly the other two parameters that we would need to tune...

%% first load some barcodes

foldL = '/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/data/Chromosome';
foldL ='/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/data/Chromosome/Yeast 280721/Raw kymos';
foldL = '/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/data/BAC/3 BACs AFF1 gene/Raw kymographs';
foldL= '/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/Fragment för matchning - Data till Erik/Mapping_New E.coli/final data/20220225_Sample-dEC-st11_110nmPERpx_0.198nmPERbp/kymos';
foldL = '/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/Fragment för matchning - Data till Erik/Mapping_New E.coli/collected';
% load single frame bars
import Core.load_single;
[barcodeGen,passingKymos,kymoStructs,kS,filelist] = load_single(foldL);

% now laod 2
fold2 ='/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/data/Chromosome/EF365 (Sample 2)/211018_Sample5.2_110nmPERpx_0.2252 nmPERbp/filter2_int_45/kymos.final_211018';
fold2 ='/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/Fragment för matchning - Data till Erik/211101_Sample 5.2_110nmPERpx_0.202nmPERbp';
[barcodeGen2,passingKymos2,kymoStructs2,kS2,filelist2] = load_single(fold2);

%%
minOverlapLen = 300;


% for the analysis, the total number of fitting attempts will be important,
% not the lengths of individual barcodes. For the analysis we take common
% length, and short barcode
import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 1000;
PSF_WIDTH_PIXELS = 300/110;
randLenMin = 400;
randLen2 = 400;

sF = 0.95:0.01:1.05;%

% minLen = 300; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGenGood = barcodeGen(barLens>randLenMin);
% barcodeGenGood(784) = [];
barLens2 = cellfun(@(x) sum(x.rawBitmask),barcodeGen2);
barcodeGenGood2 = barcodeGen2(barLens2>randLenMin);


%  calculates all MP and MPI
NN=1;
out='output';
tic
import Nullmodel.gen_scores_real;
[maxPCC,totLen] = gen_scores_real(barcodeGenGood,barcodeGenGood2,randLenMin,randLen2,minOverlapLen,NN,sF,out,SCAMP_LINE);
toc

%%
import Zeromodel.beta_ev_params;

casetest=1;a_fit2=[];n_fit2=[];
for k=1:NN
%     k
%     tic
%     import Zeromodel.beta_ev_params;
    N_FIXED = totLen(k)*(randLenMin+(k-1)*10-minOverlapLen+1);

    if casetest==1
        [parameters] = beta_ev_params(maxPCC{k}, minOverlapLen/2);
    else
        [parameters] = beta_ev_params(maxPCC{k}, 200, N_FIXED);
    end
    %     toc
    a_fit2(k) = parameters(1);
    n_fit2(k) = parameters(2);

%     [a_fit(k), n_fit(k)] = beta_ev_fit(maxPCC{k}, [2 1], [inf inf], [4 (RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS)], [true false]);

end
% a parameters seems to be constant, so we take the mean
aPar = mean(a_fit2);

c1 = polyfit(1:length(a_fit2),a_fit2,1) % fit a line.

par2 = (2*(randLenMin+10*(kk-1)-minOverlapLen)*(randLen2-minOverlapLen))./n_fit2

coef=aPar/minOverlapLen

figure,plot(n_fit2)
hold on
kk=1:NN;
plot(2*(randLenMin+10*(kk-1)+randLen2-2*minOverlapLen))

%% parameter 2.. recalculated based on fixed first parameter
% for k=1:NN
%     k
% %     tic
% %     import Zeromodel.beta_ev_params;
% %     N_FIXED = totLen(k)*(RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS+1);
% 
%     m = length(maxPCC{k});
%     cc = maxPCC{k};
%     x2 = aPar;
%     
%     denom = (-1 / m * sum(log(1 + betainc(cc.^2, 1/2, x2 / 2 - 1))) + log(2));
%     n2(k) = max(0, 1 ./ denom);
% end
%%

% now, x1 length is always 500, N2 length os 500+10*(k-1);
% need to remove overlap?
x = 2*(randLenMin+10*(kk-1)-minOverlapLen)*(randLen2-minOverlapLen)/par2;
% x = 2*(RAND_LENGTH_MIN+([1:NN]-1)*10-MIN_OVERLAP_PIXELS+1+RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS);
% figure,plot(n2test);hold on;plot(x)

% 
% n2test=n2;
% x = 3*(RAND_LENGTH_MIN+([1:100]-1)*10-MIN_OVERLAP_PIXELS+1)+RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS;

% x = 17*(RAND_LENGTH_MIN+([1:100]-1)*10-MIN_OVERLAP_PIXELS+1)+RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS;

% c = polyfit(x,n2,1)
% y_est = polyval(c,x);
% figure,plot(x,n2)
% hold on
% plot(x,y_est,'r--','LineWidth',2)
% 
% k=1;
% xx=0:0.001:1;
% import Zeromodel.beta_ev_pdf;
% [p] = beta_ev_pdf(xx, aPar, 1, y_est(k));
% f=figure,plot(xx,p)
% hold on
% histogram(maxPCC{k},'Normalization','pdf')


xx=0:0.001:1;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,aPar, 1, x);
f=figure,plot(xx,p)
hold on
histogram(maxPCC{k},'Normalization','pdf')
% set(gca,'YScale','log')
%}
%%
end