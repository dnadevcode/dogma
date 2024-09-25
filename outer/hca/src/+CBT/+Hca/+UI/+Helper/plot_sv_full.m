function [fullTable] = plot_sv_full( resultStruct,sets,name,qr,plotonlypass,bpPerPx)
    % plot_sv
    %
    % in this case allow the first barcode, the query, to be circular 
    % (since we probably want to do the same later. Or 
    %
    %
    %
    % plot the sample of simulation
    %
    %   Args:
    %      resultStruct,sets,name,qr
    %   Returns:
    %
    %   fullTable: structure with fulltable of data
    %
    %
    
%     if nargin < 4
%         input.tableColLabels = {'Start ref','Stop Ref','Start Query', 'Stop Query','Orientation','Alignment'};
%     else
%         input.tableColLabels = {'Start Query','Stop Query','Start Ref', 'Stop Ref','Orientation','Alignment'};
%     end
    if nargin < 6
        bpPerPx = 1;
        labelstr = 'Position (px)';
        str2 = 'px';
    else
        labelstr = 'Position (kbp)';
        str2 = 'kbp';
    end
%     
    if nargin < 4
        timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    else
        timestamp = 'test';
    end

    
%     try 
%         input.dataFormat = {'%i'};
%         input.data = randStruct.matchTable; 
%         latex = latexTable(input);
% 
% 
%         % save LaTex code as file
%         file = fullfile(sets.fold, strcat(timestamp,'table.txt'));
%         fd = fopen(file,'w');
%         [nrows,ncols] = size(latex);
%         for row = 1:nrows
%             fprintf(fd,'%s\n',latex{row,:});
%         end
%         fclose(fd);
%         fprintf('\n... your LaTex code has been saved \n');
% 

    % need to figure out the case when we allow circ = 1
%     f = figure;
    
    import CBT.Hca.UI.Helper.create_full_table;
    
    numBars = length(resultStruct);
%     for idx=1:length(resultStruct)
        try
        switch sets.svList
        case 1
            sets.svType = 'Insertion';
        case 2
            sets.svType = 'Invertion';
        case 3
            sets.svType = 'Repeat';
        case 4
            sets.svType = 'Translocation';
        otherwise
            error('wrong sv type');
        end
        catch
            sets.svType = num2str(1);
        end

        
        ax = subplot(1, 4,[1 2 3]);
%         ax = gca;

    ColOrd = get(ax, 'ColorOrder');
    bar1 = (resultStruct.bar1-nanmean(resultStruct.bar1))/nanstd(resultStruct.bar1);
    bar2 = (resultStruct.bar2-nanmean(resultStruct.bar2))/nanstd(resultStruct.bar2);
    res_table = resultStruct.matchTable;

    ylmax = max(bar1);
    ylmin = min(bar1);
    yumax = max(bar2+16);
    yumin = min(bar2+16);  
    hold on;

    h(1) = plot(ax, bar2, 'Color', 'black');
    hold on;
    h(2) = plot(ax, bar1+16, 'Color', 'red');

    legendNames = cell(1,size(res_table, 1));
%     end
    fullTable = [];
    N = length(bar1);    M = length(bar2);

    % problem: both bar1 and bar2 can be circularly shifted..
    for i = 1:size(res_table, 1)
        [tempTable] = create_full_table(res_table(i,:), bar1,bar2);
%         import functions.convert_matchtable;
%         tempTable = convert_matchtable(tempTable);
%         [tempTable] = create_full_table(tempTable, bar2,bar1);
%         tempTable = convert_matchtable(tempTable);
        for j=1:size(tempTable,1)
                pX = tempTable(j, [1 1 2 2 4 4 3 3]);
                pX = pX + [-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5];
                pY = [yumin yumax yumax yumin ylmax ylmin ylmin ylmax];
                pY = pY + [0 0.5 0.5 0 0 -0.5 -0.5 0 ];
%                  pY = pY + [0 1 1 0 0 -1 -1 0 ];
                if resultStruct.pass(i) == 1
                    h(i+2) = patch(pX, pY, ColOrd(1+mod(i, 7), :), 'faceAlpha', 0.1, ...
                          'edgeAlpha', 0.3, 'edgeColor', ColOrd(1+mod(i, 7), :));
                else
                    if ~plotonlypass
                       % this does not pass the thresh 
                        h(i+2) = patch(pX, pY,ColOrd(1+mod(i, 7), :), 'faceAlpha',0, ...
                              'edgeAlpha', 0.5, 'edgeColor',ColOrd(1+mod(i, 7), :),'LineStyle','--');
                    end
                    % uint8([17 17 17])
                end
        end
            
           fullTable = [ fullTable;  tempTable];
        if plotonlypass 
            if resultStruct.pass(i) == 1
                try
                     legendNames{i} = strcat(['$C_{' num2str(i) '}$ = ' num2str(resultStruct.fragmentpcc(i),3) ...
                         ', $p_{' num2str(i) '}$ = ' char(vpa(resultStruct.pval(i),3)) ...
                         ', $l_{' num2str(i) '}$=' num2str(resultStruct.lengths(i)) 'px' ]);
                catch
                   legendNames{i} = num2str(i);
                end
            end
        else
             try
                 legendNames{i} = strcat(['$C_{' num2str(i) '}$ = ' num2str(resultStruct.fragmentpcc(i),3) ...
                     ', $p_{' num2str(i) '}$ = ' char(vpa(resultStruct.pval(i),3)) ...
                     ', $l_{' num2str(i) '}$=' num2str(bpPerPx*resultStruct.lengths(i)) str2 ]);
            catch
               legendNames{i} = num2str(i);
            end           
            
        end
    end
        
    hold off;
%     if plotonlypass 
%         legend(h(logical([0 0 resultStruct.pass])),legendNames(resultStruct.pass),'Interpreter','latex','Location', 'southoutside')
%     else
%     legend(h(3:end),legendNames,'Interpreter','latex','Location', 'southoutside')
%     end
    ylabel('Shifted barcode intensities','FontSize', 10,'Interpreter','latex');
%     xlabel('Position','FontSize', 10,'Interpreter','latex');
%     title(sets.svType,'FontSize', 10,'Interpreter','latex')  ;

    ticks = 1:50/bpPerPx:2*length(bar2);
    ticksx = floor(ticks*bpPerPx);
    ax.XTick = [ticks];
    ax.XTickLabel = [ticksx];   
    xlabel( ax,labelstr,'FontSize', 10,'Interpreter','latex')

    [tempTable,barfragq,barfragr] = create_full_table(res_table,bar1,bar2,1);
%     import functions.convert_matchtable;
%     tempTable = convert_matchtable(tempTable);
%     [tempTable,barfragq2,barfragr2] = create_full_table(tempTable, bar2,bar1);
%     tempTable = convert_matchtable(tempTable);
    legend({strcat(['$\hat C_{\rm ' qr '}=$' num2str(pcc(barfragq{1},barfragr{1}'),'%0.4f')]), name},'Interpreter','latex','Location','southoutside')

    ax = subplot(1,4,4);

    for i=1:length(barfragq)
         if resultStruct.pass(i) == 1 || ~plotonlypass  
            plot( zscore(barfragq{i})+i*5,'color',ColOrd(1+mod(i, 7), :))
            hold on 
            plot(zscore(barfragr{i})+i*5,'black')
         end
    end
    ticks = 1:50/bpPerPx:2*length(bar2);
    ticksx =  resultStruct.matchTable(3)+floor(ticks*bpPerPx);
    ax.XTick = [ticks];
    ax.XTickLabel = [ticksx];   
    xlabel( ax,labelstr,'FontSize', 10,'Interpreter','latex')

    

    title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%     %
    import CBT.Hca.Core.Comparison.pcc;

% 


%     saveas(f,fullfile(sets.fold,strcat(timestamp,name)),'epsc')
%     fullfile(sets.fold,strcat(timestamp,name))
%     end

    % also want to do this comparison on a circle!!
    %
% end
% 
%    ColOrd = get(ax, 'ColorOrder');
%     barcode1 = zscore(barcode1);
%     barcode2 = zscore(barcode2);
%     ylmax = max(barcode1);
%     ylmin = min(barcode1);
%     yumax = max(barcode2+8);
%     yumin = min(barcode2+8);  
%     
%     plot(ax, barcode1, 'Color', 'black')
%     hold on;
%     plot(ax, barcode2+8, 'Color', 'black');
%     [res_table, ~] = parse_vtrace(vitResults);
%     for i = 1:size(res_table, 1)
%         pX = res_table(i, [1 1 2 2 4 4 3 3]);
%         pX = pX + [-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5];
%         pY = [yumin yumax yumax yumin ylmax ylmin ylmin ylmax];
%         pY = pY + [0 0.5 0.5 0 0 -0.5 -0.5 0 ];
%         patch(pX, pY, ColOrd(1+mod(i, 7), :), 'faceAlpha', 0.1, ...
%               'edgeAlpha', 0.3, 'edgeColor', ColOrd(1+mod(i, 7), :));
%     end
%     hold off;
% 
