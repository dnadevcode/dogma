%%
load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S1_data_updatedReal.mat')


alphav = 0.085;
alphaN = 0.15;


TT = length(lB)
NN = length(lA);
maxPCCAll = pccs;

import Zeromodel.beta_ev_params;

casetest=1;
intV = 5; % degrees of freedom for chi2 test
a_fit2 = zeros(TT,NN);n_fit2 = zeros(TT,NN);chi2Score = zeros(TT,NN);
for j=1:TT
    for k=1:NN
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
[f] = plot_evd_figure(pccs, a_fit2, n_fit2,  8, 3, lA, lB, par1,par3, par1CI, par3CI,alphav,alphaN)
% print('FIGS/FigS1.eps','-depsc','-r500');
print('FIGS/FigS1Real.eps','-depsc','-r500');

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
xlabel('$L_B$','Interpreter','latex')
ylabel('$N_{eff}/(2\cdot L_B)$','Interpreter','latex')
title('(B) $\alpha_{N}$ estimation','Interpreter','latex')
plot(lB,alphaN*ones(1,length(lB)),'--')
xlim([lB(1) lB(end)])
%
nexttile

imagesc(lB,lA,a_fit2);
title('(C) $\nu_{eff}$ heatmap','Interpreter','latex')
ylabel('$L_A$','Interpreter','latex')
xlabel('$L_B$','Interpreter','latex')
colorbar;colormap(gray) 
nexttile
imagesc(lB,lA,n_fit2);
colorbar;colormap(gray) 
title('(D) $N_{eff}$ heatmap','Interpreter','latex')
ylabel('$L_A$','Interpreter','latex')
xlabel('$L_B$','Interpreter','latex')


nexttile

iy = 1;
j = 1;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
[p2] = beta_ev_pdf(xx, alphav*lA(j), 1, alphaN*2*lB(iy));

% f=figure,
plot(xx,p,'LineWidth',2)
hold on

histogram(pccs{j,iy},'Normalization','pdf')
title('E) EVD fit histogram','Interpreter','latex')
xlim([0.3 0.95])
plot(xx,p2,'LineWidth',2)

lgd = legend({'Fitted EVD','PCC histogram','Pre-calculated parameters EVD'},'Interpreter','latex','Location','southoutside')

print('FIGS/FigS1Real.eps','-depsc','-r500');

%%
load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S2_data_updated2Real.mat')


import Zeromodel.beta_ev_params;
a_fit2 = zeros(size(maxPCC,1),size(maxPCC,2)); % first dim lA, second dim w
n_fit2 = zeros(size(maxPCC,1),size(maxPCC,2));
for j = 1:length(w)
    for k=1:length(lA)
            [parameters] = beta_ev_params(maxPCC{k,j},w(j)/2);%(randLenMin+(k-1)*gap)/2);
            a_fit2(k,j) = parameters(1);
            n_fit2(k,j) = parameters(2);
    end
end

lenA = length(lA);
lenw = length(w);

% PLOT FIGURE P-VAL fit
par1 = mean(a_fit2,1);
par1CI =  1.98*[std(a_fit2,[],1)./sqrt(lenA)];

%
ix = 11; % this is w = 300
par3 = n_fit2(:,ix)'; % can only do par3 for fixed w or fixed len2, since it depends on both. Write this in the text!

betanu = 0.065;
betaN = 0.002;


import Nullmodel.PvalScripts.plot_evd_figure_2;
[f] = plot_evd_figure_2(maxPCC, a_fit2, n_fit2, 4, ix, lA, w, par1, par3, par1CI,betanu,betaN)
print('FIGS/FigS2Real.eps','-depsc','-r500');

% 
% f=figure('Position',[1 1 800 800])
% tiledlayout(3,2,'TileSpacing','tight')
% 
% nexttile
% hold on
% plot(w,par1./w)
% 
% errorbar(w,par1./w,par1CI./w)
% xlabel('$w$','Interpreter','latex')
% ylabel('$\nu_{eff}/w$','Interpreter','latex')
% title('(A) $\beta_{\nu}$ estimation','Interpreter','latex')
% plot(w,betanu*ones(1,length(w)),'--')
% 
% nexttile
% hold on
% % plot(barL1,par3./(2.*(max(barL1-overlapL(ix),barL2-overlapL(ix)))))
% plot(barL1,par3./((2.*((barL1-w(ix)).*((barL2-w(ix)))))))
% % plot(barL1,par3./((2.*((barL1).*((barL2))))))
% 
% xlabel('$L_A$','Interpreter','latex')
% ylabel('$N_{eff}/(2(\cdot L_B-w)(L_A-w)$','Interpreter','latex') %(2(\cdot L_B-w)(L_A-w)
% title('(B) $\beta_{N}$ estimation','Interpreter','latex')
% % 
% plot(barL1,betaN*ones(1,length(barL1)),'--')
% 
% nexttile
% % shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% % longL = RAND_LENGTH_2+(0:NN-1)*gap;
% imagesc(barL1,w,a_fit2);
% title('(C) $\nu_{eff}$ heatmap','Interpreter','latex')
% colorbar;colormap(gray) 
% xlabel('$L_A$','Interpreter','latex')
% ylabel('$w$','Interpreter','latex')
% nexttile
% % shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% % longL = RAND_LENGTH_2+(0:NN-1)*gap;
% imagesc(barL1,w,n_fit2);
% colorbar;colormap(gray) 
% title('(D) $N_{eff}$ heatmap','Interpreter','latex')
% xlabel('$L_A$','Interpreter','latex')
% ylabel('$w$','Interpreter','latex')
% 
% nexttile
% % figure
% iy = 2;
% j = 11;
% % xx=0.1:0.001:0.9;
% % import Zeromodel.beta_ev_pdf;
% % [p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
% % % f=figure,
% % plot(xx,p)
% % hold on
% % histogram(maxPCC{iy,j},'Normalization','pdf')
% % xlabel('Max local PCC','Interpreter','latex')
% % title('E) EVD fit histogram','Interpreter','latex')
% % xlim([0.2 0.8])
% % lgd = legend({'Functional EVD','PCC histogram'},'Interpreter','latex','Location','southoutside')
% 
% xx=0.1:0.001:0.9;
% import Zeromodel.beta_ev_pdf;
% [p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
% [p2] = beta_ev_pdf(xx,betanu*w(j), 1, betaN*(2.*(lA(iy)-w(j)).*(lB(iy)-w(j))));
% 
% % f=figure,
% plot(xx,p,'LineWidth',2)
% hold on
% histogram(maxPCC{iy,j},'Normalization','pdf')
% xlabel('Max local PCC','Interpreter','latex')
% title('E) EVD fit histogram','Interpreter','latex')
% xlim([0.2 0.8])
% plot(xx,p2,'LineWidth',2)
% 
% 
% lgd = legend({'Fitted EVD','PCC histogram','Pre-calculated parameters EVD'},'Interpreter','latex','Location','southoutside')
% 
% 
% 
