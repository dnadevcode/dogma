function [fig1] = plot_best_pos( fig1, comparisonStruct, numBar, sets, markers,lengthBorders,scaleF)
    % plot_best_pos - pltos three maximum coefficients

    % requires from comparisonStruct: idx, pos, lengthMatch
    
    if isempty(fig1) % if figure is not called.
        fig1 = figure;
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

    p3 = plot(fig1, (pos+posShift)/scaleF,1:size(pos,1),'ob');
    p3(1).Marker = markers(1);
    try
        p3(2).Marker = markers(2);
    end
    try
        p3(3).Marker = markers(3);
    end
%     hold on
    

%     p4.Marker = markers(1);
%         p4 = plot(pos+posShift,1:size(pos,1),'ob');

  
    if  sets.genConsensus == 1
        plot(fig1,(0:100:sum(lengthBorders))/scaleF, 0.5+repmat(numBar,length(0:100:sum(lengthBorders)),1))
    end
    
    plot(fig1,lengthBorders/scaleF,zeros(1,length(lengthBorders)),'redx')
    
    % add lines for barcode lengths
    posEnd = cell2mat(cellfun(@(x) x.pos(1)+x.lengthMatch,comparisonStruct,'UniformOutput',0)');
    for i=1:size(pos,1)
        plot(fig1,([pos(i,1) posEnd(i)]+posShift(i))/scaleF,[i i],'ob-');     
    end
    
    if scaleF == 1
        xlabel(fig1,'Best position (px)','Interpreter','latex')  
    else
        xlabel(fig1,'Best position (Mbp)','Interpreter','latex')  
    end
    
    ylabel(fig1,'Barcode nr.','Interpreter','latex')
    if size(pos,2) == 1
        legend(fig1,{'$\hat C$','Theoretical barcodes seperator'},'Location','southoutside','Interpreter','latex')
    else
        legend(fig1,{'$\hat C$','$C_2$','$C_3$','Theoretical barcodes seperator'},'Location','southoutside','Interpreter','latex')
    end
    
    ylim(fig1,[0,size(pos,1)+2])

    title(fig1,'Best position');
    
end



