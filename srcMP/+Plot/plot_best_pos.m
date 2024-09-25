function [fig1] = plot_best_pos( fig1, comparisonStruct, numBar, sets, markers,lengthBorders,scaleF)
    % plot_best_pos - pltos three maximum coefficients
    %   Args:
    %
    %   ReturnsL
    %       fig1
    
    if isempty(fig1) % if figure is not called.
        fig1 = figure;
        tiledlayout(1,4);
        ax1=nexttile([1 3]);
    end
  
    if isempty(numBar)
        numBar = length(comparisonStruct);
    end
    
    if isempty(markers)
        markers = ['o';'s';'x';'+';'d';'v'];
    end
    
    if isempty(sets)
        sets.genConsensus = 0;
    end
    
    if nargin < 7 
        scaleF = 1;
    end
    
    cumLengths = [0 lengthBorders]'; % plot all the theories on x axis
    
    % which theory
    posShift = cumLengths(cellfun(@(x) x.idx,comparisonStruct)');
    
    % position along the theory
    pos = cell2mat(cellfun(@(x) x.pos(1:min(end,3)),comparisonStruct,'UniformOutput',0)');
    ornt = cell2mat(cellfun(@(x) x.or(1:min(end,3)),comparisonStruct,'UniformOutput',0)');
    score = cell2mat(cellfun(@(x) x.maxcoef(1:min(end,3)),comparisonStruct,'UniformOutput',0)');
    lengthMatch = cell2mat(cellfun(@(x) x.lengthMatch,comparisonStruct,'UniformOutput',0)');

%     p3 = plot((pos+posShift)/scaleF,1:size(pos,1),'ob');
    p3(1).Marker = markers(1);
    try
        p3(2).Marker = markers(2);
    catch
    end
    try
        p3(3).Marker = markers(3);
    catch
    end
    hold on
    

%     p4.Marker = markers(1);
%         p4 = plot(pos+posShift,1:size(pos,1),'ob');

  
    if  sets.genConsensus == 1
        plot((0:100:sum(lengthBorders))/scaleF, 0.5+repmat(numBar,length(0:100:sum(lengthBorders)),1))
    end
    
    p = zeros(1,3);
    p(1) = plot(lengthBorders/scaleF,zeros(1,length(lengthBorders)),'redx');
    
    % add lines for barcode lengths
    posEnd = cell2mat(cellfun(@(x) x.pos(1)+x.lengthMatch,comparisonStruct,'UniformOutput',0)');
    for i=1:size(pos,1)
        if ornt(i)==1
            p(2) = plot(([pos(i,1) posEnd(i)]+posShift(i))/scaleF,[i i],'>b-');     
        else
            p(3) = plot(([pos(i,1) posEnd(i)]+posShift(i))/scaleF,[i i],'<r-');     
        end
    end
    
    % we also plot colorbar for best score!
    
    if scaleF == 1
        xlabel('Best position (px)','Interpreter','latex')  
    else
        xlabel('Best position (Mbp)','Interpreter','latex')  
    end
    
    legendList = {'Islands separator','Barcodes or=1','Barcodes or=2'};
    
    ylabel('Barcode nr.','Interpreter','latex')
%     if size(pos,2) == 1
        legend(p(find(p)), legendList(find(p)) ,'Location','southoutside','Interpreter','latex')
%     else
%         legend({'$\hat C$','$C_2$','$C_3$','Theory edge seperator'},'Location','southoutside','Interpreter','latex')
%     end
    
    ylim([0,size(pos,1)+2])
   
    grid on
    grid minor
    ax2=nexttile
    scaledScore = score(:,1)./sqrt(lengthMatch);
    
    hHM =heatmap([score(:,1) scaledScore./max(scaledScore)])
    hHM.NodeChildren(3).YDir='normal';                   % turn Y-Axis normal directionend
%     hHM.XDisplayLabels = nan(size(hHM.XDisplayData));
%     hHM.YDisplayLabels = nan(size(hHM.YDisplayData));
    xlabel('score')
%     hHM.CellLabelColor='none';
%     linkaxes([ax1 ax2])


