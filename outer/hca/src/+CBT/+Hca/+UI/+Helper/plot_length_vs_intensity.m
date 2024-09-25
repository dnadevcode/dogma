function [] = plot_length_vs_intensity(ax, barcodeGenC)


DNAExtensions = cellfun(@(x) length(x.rawBarcode),barcodeGenC) ;
intensities = cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGenC);


color = {[0, 0.65, 0.95],'cyan','none'}; % magenta, cyan, grey
sz = [40,30,30];
linewidth = [0.5, 1.5,1.5];
alpha = [0.6,0.6,0.6];
ax.XLim = [0,max(DNAExtensions)+10];
ax.YLim = [0,max(intensities)*1.1];
hold(ax,'on')

s = scatter(ax,DNAExtensions,intensities,sz(1),...
    'MarkerEdgeColor',color{1},'MarkerFaceColor',color{1},...
    'Linewidth',linewidth(1),...
    'MarkerFaceAlpha',alpha(1),'MarkerEdgeAlpha',alpha(1));

% for i = 1:length(DNAExtensions)
%     xData = DNAExtensions(i);
%     yData = intensities(i);
%     s = scatter(ax,xData,yData,sz(i),...
%         'MarkerEdgeColor',color{i},'MarkerFaceColor',color{i},...
%         'Linewidth',linewidth(i),...
%         'MarkerFaceAlpha',alpha(i),'MarkerEdgeAlpha',alpha(i));
% end

% 
% s.LineWidth = 1;
% s.MarkerEdgeColor = [0.3 0.3 0.3];
% s.MarkerFaceColor = 'none';

% legend(ax,'lambda','non-lambda','location','southoutside');
xlabel(ax,'Extension (px)')
ylabel(ax,'Intensity')
title(ax,'Length vs intensity plot')

end

