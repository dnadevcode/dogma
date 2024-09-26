function [f,g,tpr] = position_evaluation_stouffer(sortedVals,sortedValsBad, foS,posShift,posShiftBad,f,g)
    % Pairwise evaluation plot position_evaluation_stouffer
    % N = 200;
    % N2 = 100;
    % msThresh = 50;

%   passthreshC{idxRun}  these passed thresh
if nargin < 6
f=figure,g=tiledlayout(5,2,'TileSpacing','tight')
end
nexttile([1 2])

% plot(sortedVals,'o','Color',[0 0.4470 0.7410])
hold on
histogram(sortedVals)
histogram(sortedValsBad)

% plot(sortedVals,'o','Color',[0.8500 0.3250 0.0980])
% plot(foS,'^','Color',[0 0.4470 0.7410])
% plot(foS,'^','Color',[0.8500 0.3250 0.0980])
grid on
    grid minor
    lgd = legend({'$s_{ \rm Stouffer}$ ','$s_{ \rm Removed}$ '},'Interpreter','latex');
    lgd.Location = 'eastoutside';
    title('(B) Stouffer scores','Interpreter','latex')
    set(gca, 'YScale', 'log')
nexttile([1 1])

% plot(posShift,'o','Color',[0.8500 0.3250 0.0980])
histboxes1 = [0:30 10000];
% histboxes2 = [histboxes1(2:end)];

h = histcounts(posShift,histboxes1);
% figure;
bar([histboxes1(1:end-1)],h,1)

set(gca,'xtick',[1 3 10 15 30.5],'xticklabel',{1,3,10, 15,'30+'}) 
xlabel('Pairwise distance (px)','Interpreter','latex')
% h2 = histogram('BinCounts', h, 'BinEdges', histboxes1);
% xlim([0 100])

% histogram(posShift, histboxes)
% set(gca, 'YScale', 'log')

% hold on
% histogram(posShiftBad)

grid on
    grid minor
    ylabel('Log of Histogram Counts','Interpreter','latex')
% lgd2 = legend({'d'},'Interpreter','latex');
lgd2.Location = 'eastoutside';
ylim([0 max(h)])
%     ylim([0 300])
% 
% xlim([0 100])
title('(C) Distance to ground truth','Interpreter','latex')
nexttile([1 1])


msThresh = 1:2:10;
valLegend = cell(1,length(msThresh));

for i=1:length(msThresh) 
    neg = posShift > msThresh(i);
    % pos
    pos = ones(1,length(posShift));
    
    % fp = neg.*(pos); % false positives
    tp = (1-neg).*pos;
    % tn = (1-pos).*neg;
    fn = (1-pos)+pos.*neg;
    
    % fpC = cumsum(fp);
    tpC = cumsum(tp);
    % tnC = cumsum(tn);
    fnC = cumsum(fn);
    
    % ppv = tpC./(fpC+tpC);
    tpr = tpC./(tpC+fnC);
    plot(tpr,'linewidth',1.5);
    hold on
    valLegend{i} = [num2str( msThresh(i))];
end
% plot(ppv)
lgd2=legend (valLegend,'Location','southoutside','Interpreter','latex');
lgd2.Location = 'eastoutside';
title(lgd2,'Threshold (px)','Interpreter','latex')
xlabel('Overlap index (sorted by $s_{Stouffer}$)','Interpreter','latex')
title('(D) True positive rate','Interpreter','latex')
ylim([0.3 1.05])

xlim([0 200])

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
% print('FIGS/FigS5.eps','-depsc','-r300');
% 

end

