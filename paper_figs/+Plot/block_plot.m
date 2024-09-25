function [consensusToPlot1] = block_plot(sI,nr,outConsensus,consensusToPlot1)


% get best indices
% [a,idx] = max(cellfun(@(x) size(x,1),outConsensus));
idx = sI(nr);

% sort based on starting position
[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
%


if nargin < 4
    consensusToPlot1 = outConsensus{idx}(idxv,:);

end


consensusToPlot = outConsensus{idx}(idxv,:);

imagesc(consensusToPlot);colormap(gray)
xlim([1 size(consensusToPlot1,2)])
ylim([0.5 size(consensusToPlot,1)+0.5])

imagesc(outConsensus{idx}(idxv,:));colormap(gray)
axis off
set(gca,'xtick',[])

end

