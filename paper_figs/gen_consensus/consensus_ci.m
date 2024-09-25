function [f] = consensus_ci(aaBarcodeIslandBlockRep,rawBarcodeIslandBlockRep,matRep)
% consensus comparison with and without re-adjustment procedure

addlim = 0;
limB = [90 180];

f=figure('Position',[1 1 1000 600]);
tiledlayout(3,1)
nexttile
title('Amplitude adjusted','Interpreter','latex')
hold on
y = mean(aaBarcodeIslandBlockRep,'omitnan');
stdy = std(aaBarcodeIslandBlockRep,1,'omitnan');
numB = sum(~isnan(aaBarcodeIslandBlockRep));

yconf = [y+1.98*stdy./sqrt(numB) y(end:-1:1)-1.98*stdy(end:-1:1)./sqrt(numB(end:-1:1))];
yconf(isnan(yconf)) = nanmean(y);
% yconf(isnan(yconf)) = y(isnan(yconf));
%
% figure

p = fill([1:length(yconf)/2 length(yconf)/2:-1:1],yconf,'red');
p.FaceColor = [1 0.8 0.8];
p.EdgeColor = 'none';
hold on
plot(1:length(y),y)
set(gca, 'YScale', 'log')
% hold on
% plot(x,y,'r-')
if addlim
ylim(limB)
end

nexttile
title('Without amplitude adjustment','Interpreter','latex')
hold on
y = mean(rawBarcodeIslandBlockRep,'omitnan');
stdy = std(rawBarcodeIslandBlockRep,1,'omitnan');
numB = sum(~isnan(rawBarcodeIslandBlockRep));

yconf = [y+1.98*stdy./sqrt(numB) y(end:-1:1)-1.98*stdy(end:-1:1)./sqrt(numB(end:-1:1))];
yconf(isnan(yconf)) = nanmean(y);
% yconf(isnan(yconf)) = y(isnan(yconf));
%
% figure

p = fill([1:length(yconf)/2 length(yconf)/2:-1:1],yconf,'red');
p.FaceColor = [1 0.8 0.8];
p.EdgeColor = 'none';
hold on
plot(1:length(y),y)
set(gca, 'YScale', 'log')
if addlim
ylim(limB)
end
% rawBarcodeIslandBlockRep
if nargin >=3
    nexttile
    title('Amplitude adjusted and re-aligned','Interpreter','latex')
    hold on
    y = mean(matRep,'omitnan');
    stdy = std(matRep,1,'omitnan');
    numB = sum(~isnan(matRep));

    yconf = [y+1.98*stdy./sqrt(numB) y(end:-1:1)-1.98*stdy(end:-1:1)./sqrt(numB(end:-1:1))];
    yconf(isnan(yconf)) = nanmean(y);
    % yconf(isnan(yconf)) = y(isnan(yconf));
    %
    % figure

    p = fill([1:length(yconf)/2 length(yconf)/2:-1:1],yconf,'red');
    p.FaceColor = [1 0.8 0.8];
    p.EdgeColor = 'none';
    hold on
    plot(1:length(y),y)
    set(gca, 'YScale', 'log')
    xlabel('Pixel position','Interpreter','latex')
    legend({'95% CI','Consensus'},'location','southoutside','Interpreter','latex')
    if addlim

    ylim(limB)
    end
end

