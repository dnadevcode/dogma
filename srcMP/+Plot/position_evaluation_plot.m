function [ppv,tpr] = position_evaluation_plot(passLocalThresh,passGlobalThresh,sortedVals,foS,posShift, N)

% N = 200;
% N2 = 100;

%   passthreshC{idxRun}  these passed thresh
figure,tiledlayout(3,1,'TileSpacing','tight')
nexttile
plot(sortedVals(1:N),'o','Color',[0 0.4470 0.7410])
hold on
plot(find(passLocalThresh),sortedVals(passLocalThresh),'o','Color',[0.8500 0.3250 0.0980])
plot(foS(1:N),'^','Color',[0 0.4470 0.7410])
plot(find(passGlobalThresh),foS(passGlobalThresh),'^','Color',[0.8500 0.3250 0.0980])
grid on
    grid minor
    lgd = legend({'$$C $$ local,bad ','$$C $$ local, good','$$CC$$ global, bad','$$CC$$  global, good'},'Interpreter','latex')
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

msThresh = 30;
neg = posShift > msThresh;
% pos
pos = passLocalThresh'.*passGlobalThresh';

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

