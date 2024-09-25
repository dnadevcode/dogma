function [] = plot_bloc_rep(outConsensus,tsVisual)

      tsBlock = uitabgroup('Parent',tsVisual.block,'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);

      hblock = cell(1,length(outConsensus));
      tb = cell(1,length(outConsensus));
    for idx=1:length(outConsensus)
        hblock{idx} = uitab(tsBlock, 'title', strcat(['Block rep_' num2str(idx)]));
%         tiledlayout(hPanelPlot,2,1,'TileSpacing','tight','Padding','tight');
%         nexttile([1 1])
        tb{idx} = tiledlayout(hblock{idx} ,2,1,'TileSpacing','tight','Padding','tight');
      
%     figure,
%         tiledlayout(2,1);
        nexttile(tb{idx})
        hold on
        title(['Block representation (', num2str(idx),')'], 'Interpreter','latex')

        % sort based on starting position
        [pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus{idx}(x,:)),1,'first') +(find(~isnan(outConsensus{idx}(x,:)),1,'last')-find(~isnan(outConsensus{idx}(x,:)),1,'first'))/2,1:size(outConsensus{idx},1)));
        % 
        consensusToPlot1 = outConsensus{idx}(idxv,:);
        
        imagesc(consensusToPlot1);colormap(gray)
        % xlim([1 size(consensusToPlot1,2)])
        ylim([0 size(consensusToPlot1,1)+1])
        
        imagesc(outConsensus{idx}(idxv,:));colormap(gray)
        axis off
        
        set(gca,'xtick',[])
                    
        %             xlim([1 maxConsensus])
    end
            %


end

