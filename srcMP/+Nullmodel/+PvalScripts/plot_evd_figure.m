function [f] = plot_evd_figure(pccs, a_fit2, n_fit2, iy, j, lA, lB, par1, par3, par1CI, par3CI,alphav,alphaN)

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

% iy = 14;
% j = 10;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(j,iy), 1, n_fit2(j,iy));
[p2] = beta_ev_pdf(xx, alphav*lA(iy), 1, alphaN*2*lB(j));

% f=figure,
plot(xx,p,'LineWidth',2)
hold on

histogram(pccs{j,iy},'Normalization','pdf')
title(['(E) EVD fit histogram,', ' $L_A=$ ',num2str(lA(iy)),' $L_B=$ ',num2str(lB(j))  ],'Interpreter','latex')
xlim([min(pccs{j,iy})-0.1 max(pccs{j,iy})+0.1])
plot(xx,p2,'LineWidth',2)

lgd = legend({'Fitted EVD','PCC histogram','Pre-calculated parameters EVD'},'Interpreter','latex','Location','southoutside')


end

