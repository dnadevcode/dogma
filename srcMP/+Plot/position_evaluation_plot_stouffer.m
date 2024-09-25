function [ppv,tpr] = position_evaluation_plot_stouffer(sortedVals,foS,posShift,posShiftBad)

% N = 200;
% N2 = 100;
msThresh = 30;

%   passthreshC{idxRun}  these passed thresh
figure,tiledlayout(3,1,'TileSpacing','tight')
nexttile
% plot(sortedVals,'o','Color',[0 0.4470 0.7410])
hold on
plot(sortedVals,'o','Color',[0.8500 0.3250 0.0980])
% plot(foS,'^','Color',[0 0.4470 0.7410])
% plot(foS,'^','Color',[0.8500 0.3250 0.0980])
grid on
    grid minor
    lgd = legend({'Stouffer score plot'},'Interpreter','latex')
    lgd.Location = 'eastoutside';
    title('Local and global scores')
nexttile
plot(posShift,'o','Color',[0.8500 0.3250 0.0980])
grid on
    grid minor
lgd2 = legend({'Pairwise distance'});
lgd2.Location = 'eastoutside';
title('Pairwise distance plot')
nexttile

neg = posShift > msThresh;
% pos
pos = ones(1,length(posShift));

fp = neg.*(pos); % false positives
tp = (1-neg).*pos;
tn = (1-pos).*neg;
fn = (1-pos)+pos.*neg;

fpC = cumsum(fp);
tpC = cumsum(tp);
tnC = cumsum(tn);
fnC = cumsum(fn);

ppv = tpC./(fpC+tpC);
tpr = tpC./(tpC+fnC);
plot(tpr);
hold on
plot(ppv)
lgd2=legend ({'Recall','Precision'},'Location','southoutside');
lgd2.Location = 'eastoutside';

title('Precision and recall')



% hold on
% plot(foS(1:N2),'.')
% nexttile
% plot(posShift,'x');hold on
% plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
% legend({'Distance between estimated GT and pair overlap','Threshold for good scores'},'Location','southoutside')
% nexttile
% plot(posShift(1:100));hold on
% plot([length(barsPassThreshMP) length(barsPassThreshMP)], [0 max(posShift)],'red')
% legend({'Distance between estimated GT and pair overlap','Threshold for good scores'},'Location','southoutside')
print('FIGS/FigS5.eps','-depsc','-r300');


end

