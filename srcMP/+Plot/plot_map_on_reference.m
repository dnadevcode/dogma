function [fig1] = plot_map_on_reference( fig1, comparisonStruct, markers,lengthBorders,scaleF)
    % plot_map_on_reference - plots map on reference
    %   Args:
    %
    %   ReturnsL
    %       fig1
    
    if isempty(fig1) % if figure is not called.
%         fig1 = figure;
%         tiledlayout(1,4);
%         ax1=nexttile([1 3]);
    end
  
    numBar = length(comparisonStruct);
    
    if isempty(markers)
        markers = ['o';'s';'x';'+';'d';'v'];
    end
    
%     if isempty(sets)
%         sets.genConsensus = 0;
%     end
    
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

  
%     if  sets.genConsensus == 1
%         plot((0:100:sum(lengthBorders))/scaleF, 0.5+repmat(numBar,length(0:100:sum(lengthBorders)),1))
%     end
    
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
    
    legendList = {'Edge of theory barcode','Barcodes or=1','Barcodes or=2'};
    
    ylabel('Barcode nr.','Interpreter','latex')
    legend(p(find(p)), legendList(find(p)) ,'Location','southoutside','Interpreter','latex')

    
    ylim([0,size(pos,1)+2])
   
    grid on
    grid minor


