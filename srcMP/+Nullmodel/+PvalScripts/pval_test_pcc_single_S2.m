%% p-val test single with the same parameters as figure Fig4.

outmat = 'S2_data_equal_prob_bp.mat';

precalc_pval_test_pcc_single_S2(outmat)
% save('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S2_data_updated2.mat','maxPCC','osOutput','lA','w');

%%

load('/export/scratch/albertas/data_temp/bargrouping/PAPER_DATA/null/S2_data_equal_prob_bp.mat')


import Zeromodel.beta_ev_params;

casetest=1;
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

% par3 = mean(n_fit2,2);
% par3CI =  1.98*[std(n_fit2,[],2)./sqrt(lenw)];

% Similar plot for local pcc
% barL2 = lA;
% barL1 = lA;

%
ix = 11; % this is w = 300
par3 = n_fit2(:,ix)';% can only do par3 for fixed w or fixed len2, since it depends on both. Write this in the text!
overlapL = w(ix);
% minL = randLenMin+gap*([1:TT]-1);

betanu = 0.085;
betaN = 0.004;

import Nullmodel.PvalScripts.plot_evd_figure_2;
[f] = plot_evd_figure_2(maxPCC, a_fit2, n_fit2, 4, ix, lA, w, par1, par3, par1CI,betanu,betaN)

% 
% f=figure('Position',[1 1 800 800])
% tiledlayout(3,2,'TileSpacing','tight')
% 
% nexttile
% hold on
% % plot(w,par1./w)
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
% plot(barL1,par3./((2.*((lA-w(ix)).*((lA-w(ix)))))))
% % plot(barL1,par3./((2.*((barL1).*((barL2))))))
% 
% xlabel('$L_A$','Interpreter','latex')
% ylabel('$N_{eff}/(2(\cdot L_B-w)(L_A-w)$','Interpreter','latex') %(2(\cdot L_B-w)(L_A-w)
% title('(B) $\beta_{N}$ estimation','Interpreter','latex')
% % 
% plot(lA,betaN*ones(1,length(lA)),'--')
% 
% nexttile
% % shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% % longL = RAND_LENGTH_2+(0:NN-1)*gap;
% imagesc(w,lA,a_fit2);
% title('(C) $\nu_{eff}$ heatmap','Interpreter','latex')
% colorbar;colormap(gray) 
% xlabel('$w$','Interpreter','latex')
% ylabel('$L_A$','Interpreter','latex')
% nexttile
% % shortL = RAND_LENGTH_MIN+(0:TT-1)*gap;
% % longL = RAND_LENGTH_2+(0:NN-1)*gap;
% imagesc(lA,w,n_fit2');
% colorbar;colormap(gray) 
% title('(D) $N_{eff}$ heatmap','Interpreter','latex')
% xlabel('$L_A$','Interpreter','latex')
% ylabel('$w$','Interpreter','latex')
% 
% nexttile
% % figure
% iy = 4;
% j = 11;
% xx=0.1:0.001:0.9;
% import Zeromodel.beta_ev_pdf;
% [p] = beta_ev_pdf(xx,a_fit2(iy,j), 1, n_fit2(iy,j));
% [p2] = beta_ev_pdf(xx,betanu*w(j), 1, betaN*(2.*(lA(iy)-w(j)).*(lA(iy)-w(j))));
% 
% % f=figure,
% plot(xx,p,'LineWidth',2)
% hold on
% histogram(maxPCC{iy,j},'Normalization','pdf')
% xlabel('Max local PCC','Interpreter','latex')
% title(['E) EVD fit histogram,', ' $L_A=$ ',num2str(lA(iy)),', $w=$ ',num2str(w(j)) ],'Interpreter','latex')
% xlim([0.4 0.8])
% plot(xx,p2,'LineWidth',2)
% 
% % lgd = legend({'Functional EVD','PCC histogram'},'Interpreter','latex','Location','southoutside')
% lgd = legend({'Fitted EVD','PCC histogram','Pre-calculated parameters EVD'},'Interpreter','latex','Location','southoutside')
% %

print('FIGS/FigS2.eps','-depsc','-r500');


%%
alphaNu = a_fit2(end,end);
pthresh = 0.01;

    % pval MP
    alphaN = par2(end,end); %should be rather insensitive to specific value, check figure S5-S6;
    import Zeromodel.beta_ev_cdf; % correct form?
%     pvalfun = @(x,l1,l2,nuF,w) 1-beta_ev_cdf(x,nuF*w,1,nuF*2*(max(l1,l2)-w+1),0);
    pvalfun = @(x) 1-beta_ev_cdf(x,alphaNu*minOverlap,1,alphaN*2*(curLen-minOverlap+1)*(curLen-minOverlap+1),0);

    pvalLocal= pvalfun(maxPCC{end,end} )


% par1 = zeros(1,length(w));
% par2 = zeros(1,length(w));
% 
% for j = 1:length(w)
%     
%     import Zeromodel.beta_ev_params;
%     [parameters] = beta_ev_params(maxPCC{j},w(j)/2);%(randLenMin+(k-1)*gap)/2);
%     
%     par1(j) = parameters(1)/w(j);
%     par2(j) = parameters(2)./((2.*((curLen-w(j)).*((curLen-w(j))))));
% 
% end
% %% Test plot

% 
% % %     
% % [sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver, pvalCombined, sortedValsBad, sortedIdsBad] = ...
% %     calculate_sorted_pvals(oS,minOverlap, alphaNu, pthresh); %localCCStruct
% 
% 
%     stoufMP = norminv(1-pvalLocal);
% %     stoufLeft = norminv(1-pvalLeftOver);
% 
%     figure,histogram(stoufMP);hold on
% %     histogram(stoufLeft)
