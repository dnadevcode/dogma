%% p-val test single with the same parameters as figure Fig4.

% pre-generate data for figure.
precalc_pval_test_pcc_single_simple_S1

%%
% where to draw lines in plots with selected values
alphav = 0.09;
alphaN = 0.42;

% preload calculated data
load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S1_data_equal_prob_bp.mat')
lenB = length(lB);
lenA = length(lA);
maxPCCAll = pccs;

import Zeromodel.beta_ev_params;
a_fit2 = zeros(lenB,lenA);n_fit2 = zeros(lenB,lenA);chi2Score = zeros(lenB,lenA);
for j=1:lenB
    for k=1:lenA
        [parameters] = beta_ev_params(maxPCCAll{j,k}, lA(k)/2);
        a_fit2(j,k) = parameters(1);
        n_fit2(j,k) = parameters(2);
    end
end

% PLOT FIGURE P-VAL fit
par1 = mean(a_fit2,1);
par1CI =  1.98*[std(a_fit2,[],1)./sqrt(lenB)];

par3 = mean(n_fit2,2);
par3CI =  1.98*[std(n_fit2,[],2)./sqrt(lenA)];

import Nullmodel.PvalScripts.plot_evd_figure;
[f] = plot_evd_figure(pccs, a_fit2, n_fit2,  14, 10, lA, lB, par1,par3, par1CI, par3CI,alphav,alphaN)
print('FIGS/FigS1.eps','-depsc','-r500');

% plot

%%
% plot
f=figure('Position',[1 1 800 800])
tiledlayout(3,2,'TileSpacing','tight')

nexttile
hold on
% plot(lA,par1./lA)
errorbar(lA,par1./lA,par1CI./lA)

xlabel('$L_A$','Interpreter','latex')
ylabel('$\nu_{eff}/L_A$','Interpreter','latex')
title('(A) $\alpha_{\nu}$ estimation','Interpreter','latex')
plot(lA,alphav*ones(1,length(lA)),'--')
xlim([lA(1) lA(end)])

nexttile
hold on


errorbar(lB,par3'./(2*lB),par3CI'./(2*lB))
% plot(lB,par3./(2*lB))
xlabel('$L_B$','Interpreter','latex')
ylabel('$N_{eff}/(2\cdot L_B)$','Interpreter','latex')
title('(B) $\alpha_{N}$ estimation','Interpreter','latex')
plot(lB,alphaN*ones(1,length(lB)),'--')
xlim([lB(1) lB(end)])
%
nexttile

imagesc(lA,lB,a_fit2);
title('(C) $\nu_{eff}$ heatmap','Interpreter','latex')
ylabel('$L_B$','Interpreter','latex')
xlabel('$L_A$','Interpreter','latex')
colorbar;colormap(gray) 
nexttile
imagesc(lB,lA,n_fit2');
colorbar;colormap(gray) 
title('(D) $N_{eff}$ heatmap','Interpreter','latex')
ylabel('$L_A$','Interpreter','latex')
xlabel('$L_B$','Interpreter','latex')


nexttile

iy = 14;
j = 10;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
[p2] = beta_ev_pdf(xx, alphav*lA(iy), 1, alphaN*2*lB(j));

% f=figure,
plot(xx,p,'LineWidth',2)
hold on

histogram(pccs{j,iy},'Normalization','pdf')
title(['(E) EVD fit histogram,', ' $L_A=$ ',num2str(lA(iy)),' $L_B=$ ',num2str(lB(j))  ],'Interpreter','latex')
xlim([0.35 0.6])
plot(xx,p2,'LineWidth',2)

lgd = legend({'Fitted EVD','PCC histogram','Pre-calculated parameters EVD'},'Interpreter','latex','Location','southoutside')

