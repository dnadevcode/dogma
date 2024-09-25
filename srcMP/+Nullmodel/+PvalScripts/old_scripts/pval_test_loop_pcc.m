% test: generate p-values for bargroups
% SCAMP_LINE = '~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
% delete(gcp('nocreate'))
% numWorkers = 32;
% 
% parpool('local',numWorkers)
% This will depend on these parameters: overlap length, number of attempts, re-scaling factors, psf.
% possibly re-scaling factors and psf have minor relevance, then it would
% be mainly the other two parameters that we would need to tune.
%%
rng("default")
% for the analysis, the total number of fitting attempts will be important,
% not the lengths of individual barcodes. For the analysis we take common
% length, and short barcode
import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 100;
PSF_WIDTH_PIXELS = 300/110;
RAND_LENGTH_MIN = 600;
RAND_LENGTH_2 = 4000;
minOverlapLen = RAND_LENGTH_MIN;
% overlap length in general will be fixed for the analysis.

sF = 1;%
gap = 100;
% sF = 0.95:0.01:1.05;

TT = 30;
%  calculates all MP and MPI
NN = 30;

nufitEVD  =cell(1,TT);
maxPCCAll = cell(1,TT);
allPCC =  cell(1,TT);
for kk=1:TT
    
    out ='output';
    tic
    import Nullmodel.gen_scores_pcc;
    [maxPCCAll{kk},allPCC{kk}] = gen_scores_pcc(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN+(kk-1)*gap,RAND_LENGTH_2,NN,sF,out, gap);
    toc

end

% maxPCC = maxPCCAll{1};
%%
import Zeromodel.beta_ev_params;

casetest=1;
intV = 5; % degrees of freedom for chi2 test
a_fit2 = zeros(TT,NN);n_fit2 = zeros(TT,NN);chi2Score = zeros(TT,NN);
for j=1:TT
    for k=1:NN
        if casetest==1
            [parameters] = beta_ev_params(maxPCCAll{j}{k}, (RAND_LENGTH_MIN+(j-1)*gap)/2);
        else
            N_FIXED = totLen(k)*(randLenMin+(k-1)*10-minOverlapLen+1);
            [parameters] = beta_ev_params(maxPCCAll{j}{k}, 200, N_FIXED);
        end
        a_fit2(j,k) = parameters(1);
        n_fit2(j,k) = parameters(2);
            %     [a_fit(k), n_fit(k)] = beta_ev_fit(maxPCC{k}, [2 1], [inf inf], [4 (RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS)], [true false]);
    
            % Run chi2 test to see if we should keep this pair of scores
%         [histAll,vals] = histcounts(maxPCCAll{j}{k},30,'Normalization','count');
%         bincenters = (vals(2:end)+vals(1:end-1))/2;
%         cc = bincenters;
%         [pdfPlot] = beta_ev_pdf(cc,a_fit2(j,k), 1, n_fit2(j,k));
%     % 
%     %     pdfPlot = arrayfun(@(x) pdfF(x, [parameters(i,1)  parameters(i,2)]),cc);
%     %     % 
%     %     % % implement chi^2 test for the goodness-of-fit
%         import Zeromodel.chi2_test;
%         [chi2Score(j,k)] = chi2_test(histAll, pdfPlot,[1 length(histAll)],intV); % make a histogram of rand scores?
    %     pval = chi2pdf(chi2Score(j,k),intV-2 );%/length(intVals); 
    % 
    % %     p = 1-chi2cdf(chi2Score(j,k),df); % pval
    %     if pval > 0.1
    %         chi2Score(j,k) = nan;
    %     end
    end
end


%% PLOT FIGURE P-VAL fit

par1 = zeros(1,TT);
par2 = zeros(1,TT);

% choose 1
for j = 1:TT;

%     aPar = a_fit2(j,:)/minOverlapLen;
    
    kk=1:NN;
    nPar = n_fit2(j,:);
    parEq = (2*(RAND_LENGTH_2+gap*(kk-1)));
    
    % fnpar = polyfit(parEq,nPar,1) % fit a line.    
    mdl = fitlm(parEq,nPar,'intercept',false);
    fnpar(1) = mdl.Coefficients.Estimate(1);
    fnpar(2) = 0;

    par2(j) = [mdl.Coefficients.Estimate(1)];
    par1(j) = [mean(a_fit2(j,:))];
end

par3 =  zeros(1,NN);
for j=1:NN
    nPar = n_fit2(:,j);
    par3(j) = mean(nPar);
end


minL = RAND_LENGTH_MIN+gap*([1:TT]-1);
longL = RAND_LENGTH_2+gap*([1:NN]-1);

alphav = 0.14;
alphaN = 0.2;

f=figure('Position',[1 1 800 800])
tiledlayout(3,2,'TileSpacing','tight')

nexttile
hold on
plot(minL,par1./minL)
xlabel('$L_A$','Interpreter','latex')
ylabel('$\nu_{eff}/L_A$','Interpreter','latex')
title('(A) $\alpha_{\nu}$ estimation','Interpreter','latex')
plot(minL,alphav*ones(1,length(minL)),'--')
nexttile
hold on
plot(longL,par3./(2*longL))
xlabel('$L_B$','Interpreter','latex')
ylabel('$N_{eff}/(2\cdot L_B)$','Interpreter','latex')
title('(B) $\alpha_{N}$ estimation','Interpreter','latex')
plot(longL,alphaN*ones(1,length(longL)),'--')

nexttile
shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
longL = RAND_LENGTH_2+(0:NN-1)*gap;
imagesc(longL,shortL,a_fit2);
title('(C) $\nu_{eff}$ heatmap','Interpreter','latex')
ylabel('$L_A$','Interpreter','latex')
xlabel('$L_B$','Interpreter','latex')
colorbar;colormap(gray) 
nexttile
shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
longL = RAND_LENGTH_2+(0:NN-1)*gap;
imagesc(longL,shortL,n_fit2);
colorbar;colormap(gray) 
title('(D) $N_{eff}$ heatmap','Interpreter','latex')
ylabel('$L_A$','Interpreter','latex')
xlabel('$L_B$','Interpreter','latex')

nexttile

iy = 30;
j = 1;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
% f=figure,
plot(xx,p)
hold on
histogram(maxPCCAll{j}{iy},'Normalization','pdf')
title('E) EVD fit histogram','Interpreter','latex')
xlim([0.2 0.6])
lgd = legend({'Functional EVD','PCC histogram'},'Interpreter','latex','Location','southoutside')

print('FIGS/FigS1.eps','-depsc','-r300');

%%

% parEq = parEq*mean(par2);
% plot(min(x):max(x),mdl.Coefficients.Estimate(1)*(min(x):max(x)))
% lgd=legend({['$R^2$ =',num2str(mdl.Rsquared.Ordinary)],['f(x) = ',num2str(mdl.Coefficients.Estimate(1))]},'Interpreter','latex')

% f=figure,
% tiledlayout(1,2)
% nexttile
% shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% longL = RAND_LENGTH_2+(0:NN-1)*gap;
% imagesc(longL,shortL,a_fit2);
% title('$\nu_{eff}$ heatmap','Interpreter','latex')
% nexttile
% shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% longL = RAND_LENGTH_2+(0:NN-1)*gap;
% imagesc(longL,shortL,n_fit2);
% title('$N_{eff}$ heatmap','Interpreter','latex')
% ylabel('Length short')
% xlabel('Length long')

% 
% %
% saveas(f,fullfile(foldsave,'figTH.png'))
% Plot 3x2
f=figure,
tiledlayout(2,2,'TileSpacing','tight')

nexttile
plot(parEq,nPar,'x')
hold on
title('A) EVD fit $N_{eff}$','Interpreter','latex')

plot(parEq, fnpar(1)*parEq+fnpar(2))
xlabel('$v_{eff}\cdot N$','Interpreter','latex')
ylabel('$N_{eff}$','Interpreter','latex')
lgd = legend({'$N_{eff}$',['Line f(x) = ',num2str(fnpar(2)),' + ', num2str(fnpar(1)), 'x']},'Interpreter','latex');
lgd.Location = 'southoutside';
nexttile
c1 = polyfit(1:length(aPar),aPar,1) % fit a line.
coef = aPar/minOverlapLen;
plot(parEq,aPar)
hold on
plot(parEq,c1(1)*(1:length(aPar)) +c1(2))
lgd2 = legend({'$C_{\nu_{eff}}$',['Line f(x) = ',num2str(c1(2)),' + ', num2str(c1(1)), 'x']},'Interpreter','latex');
lgd2.Location = 'southoutside';
title('B) EVD fit  $\nu_{eff}$','Interpreter','latex')


nexttile

iy=1;
ix=2;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
% f=figure,
plot(xx,p)
hold on
histogram(maxPCCAll{j}{iy},'Normalization','pdf')
title('C) EVD fit histogram','Interpreter','latex')

% saveas(f,fullfile(foldsave,'figT1.png'))
%%
%%
figure
iy=6;
% % ix=2;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(iy,1), 1, n_fit2(iy,1));
% f=figure,
plot(xx,p)
hold on
histogram(maxPCC{iy},'Normalization','pdf')
title('EVD fit histogram','Interpreter','latex')


%%
% 
% %% now fit p-val dist. on these. just some tests
% import Zeromodel.beta_ev_fit;
% tic
% [a_fit, n_fit] = beta_ev_fit(maxPCC{1}, [4 1], [inf inf], [4 1], false(1, 2));
% toc
% k=10;
% 
% totLenRescaledCur = totLen(k);
% import Zeromodel.beta_ev_params;
% 
% N_FIXED = totLenRescaledCur*(RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS+1);
% [parameters] = beta_ev_params(maxPCC{k}, 200, N_FIXED);
%     
% [parameters] = beta_ev_params(maxPCC{k}, 200);
% parameters(1)
% % plot fit vs beta ev pdf
% k=5;
% xx=0:0.001:1;
% import Zeromodel.beta_ev_pdf;
% [p] = beta_ev_pdf(xx,  0.1572*MIN_OVERLAP_PIXELS, 1, 2*(500+10*(k-1)+500-2*MIN_OVERLAP_PIXELS));
% figure,plot(xx,p)
% hold on
% histogram(maxPCC{k},'Normalization','pdf')

% now calculate params actually
%      [parameters] = beta_ev_params(maxPCC{k}, MIN_OVERLAP_PIXELS/2);
import Zeromodel.beta_ev_params;

casetest=1;a_fit2=[];n_fit2=[];
for k=1:NN
    [parameters] = beta_ev_params(maxPCC{k}, RAND_LENGTH_MIN/2);
    %     toc
    a_fit2(k) = parameters(1);
    n_fit2(k) = parameters(2);

%     [a_fit(k), n_fit(k)] = beta_ev_fit(maxPCC{k}, [2 1], [inf inf], [4 (RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS)], [true false]);

end
% a parameters seems to be constant, so we take the mean
% aPar = mean(a_fit2);

kk=1:NN;
lenoverlap = RAND_LENGTH_MIN; %RAND_LENGTH_MIN+10*(kk-1);

mean(a_fit2./lenoverlap)

c1 = polyfit(1:length(a_fit2),a_fit2,1); % fit a line.


% coef=aPar/MIN_OVERLAP_PIXELS

figure,plot(n_fit2)
hold on
kk=1:NN;
plot(2*(RAND_LENGTH_2+10*(kk-1)-RAND_LENGTH_MIN))



% plot(2*(RAND_LENGTH_MIN+10*(kk-1)+RAND_LENGTH_2-2*MIN_OVERLAP_PIXELS))
% 
% par2 = (2*(RAND_LENGTH_MIN+10*(kk-1)-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_2-MIN_OVERLAP_PIXELS))./n_fit2
% mean(par2)

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
x = 2*(RAND_LENGTH_MIN+10*(kk-1)-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_2-MIN_OVERLAP_PIXELS)/par2;
% figure,plot(n_fit2);hold on;plot(x)

% 
% n2test=n2;
% x = 3*(RAND_LENGTH_MIN+([1:100]-1)*10-MIN_OVERLAP_PIXELS+1)+RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS;

% x = 17*(RAND_LENGTH_MIN+([1:100]-1)*10-MIN_OVERLAP_PIXELS+1)+RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS;
% 
% c = polyfit(x,n_fit2,1)
% y_est = polyval(c,x);
% figure,plot(x,n_fit2)
% hold on
% plot(x,y_est,'r--','LineWidth',2)
%
k=1;
xx=0:0.001:1;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx, a_fit2(k), 1, n_fit2(k));%y_est(k)
f=figure,plot(xx,p)
hold on
histogram(maxPCC{k},'Normalization','pdf')
% set(gca,'YScale','log')

intrestingPCC = maxPCC{k}(maxPCC{k}>0.7);
%%
%     N_FIXED = (RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS+1)*(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS+1)*11;
% 
% tic
%  [parameters] = beta_ev_params(maxPCC{1}, 100, N_FIXED);
% parameters(1)
% toc
% % parameters(2)
% tic
%     [a1, n1] = beta_ev_fit(maxPCC{1}, [2 1], [inf inf], [40 N_FIXED], [false true])
% toc
%     [a1, n1] = beta_ev_fit(maxPCC{1}, [2 1], [inf inf], [40 N_FIXED], [false false])

    
% figure,plot([RAND_LENGTH_MIN:10:RAND_LENGTH_MIN+(NN-1)*10]-MIN_OVERLAP_PIXELS+(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS)*length(sF),n_fit2)
% 
% cellfun(@(x) nanmean(x),maxPCC)
% 
% cell2mat(maxPCC)
% 
% 
% pccMat = cell2mat(maxPCC');


% import TempStuff.gen_random_fragments
% [barcodes, refBarcode] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD);

%% TEST: pval for real data. Here we take barcodes from unrelated bacteria, and calculate the same.

%%