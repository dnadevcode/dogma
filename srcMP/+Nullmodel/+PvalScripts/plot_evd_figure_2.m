function [f] = plot_evd_figure_2(maxPCC, a_fit2, n_fit2, iy, j, lA, w, par1, par3, par1CI,betanu,betaN)

% function [f] = plot_evd_figure(pccs, a_fit2, n_fit2, iy, j, lA, lB, par1, par3, par1CI, par3CI,alphav,alphaN)

f=figure('Position',[1 1 800 800])
tiledlayout(3,2,'TileSpacing','tight')

nexttile
hold on
% plot(w,par1./w)

errorbar(w,par1./w,par1CI./w)
xlabel('$w$','Interpreter','latex')
ylabel('$\nu_{eff}/w$','Interpreter','latex')
title('(A) $\beta_{\nu}$ estimation','Interpreter','latex')
plot(w,betanu*ones(1,length(w)),'--')

nexttile
hold on
% plot(barL1,par3./(2.*(max(barL1-overlapL(ix),barL2-overlapL(ix)))))
plot(lA,par3./((2.*((lA-w(j)).*((lA-w(j)))))))
% plot(barL1,par3./((2.*((barL1).*((barL2))))))

xlabel('$L_A$','Interpreter','latex')
ylabel('$N_{eff}/(2(\cdot L_B-w)(L_A-w)$','Interpreter','latex') %(2(\cdot L_B-w)(L_A-w)
title('(B) $\beta_{N}$ estimation','Interpreter','latex')
% 
plot(lA,betaN*ones(1,length(lA)),'--')

nexttile
% shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% longL = RAND_LENGTH_2+(0:NN-1)*gap;
imagesc(w,lA,a_fit2);
title('(C) $\nu_{eff}$ heatmap','Interpreter','latex')
colorbar;colormap(gray) 
xlabel('$w$','Interpreter','latex')
ylabel('$L_A$','Interpreter','latex')
nexttile
% shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% longL = RAND_LENGTH_2+(0:NN-1)*gap;
imagesc(lA,w,n_fit2');
colorbar;colormap(gray) 
title('(D) $N_{eff}$ heatmap','Interpreter','latex')
xlabel('$L_A$','Interpreter','latex')
ylabel('$w$','Interpreter','latex')

nexttile
% figure
% iy = 4;
% j = 11;
xx=0.1:0.001:0.9;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx,a_fit2(iy,j), 1, n_fit2(iy,j));
[p2] = beta_ev_pdf(xx,betanu*w(j), 1, betaN*(2.*(lA(iy)-w(j)).*(lA(iy)-w(j))));

% f=figure,
plot(xx,p,'LineWidth',2)
hold on
histogram(maxPCC{iy,j},'Normalization','pdf')
xlabel('Max local PCC','Interpreter','latex')
title(['E) EVD fit histogram,', ' $L_A=$ ',num2str(lA(iy)),', $w=$ ',num2str(w(j)) ],'Interpreter','latex')
% xlim([0.4 0.8])
xlim([min(maxPCC{iy,j})-0.1 max(maxPCC{iy,j})+0.1])

plot(xx,p2,'LineWidth',2)

% lgd = legend({'Functional EVD','PCC histogram'},'Interpreter','latex','Location','southoutside')
lgd = legend({'Fitted EVD','PCC histogram','Pre-calculated parameters EVD'},'Interpreter','latex','Location','southoutside');
%


end

