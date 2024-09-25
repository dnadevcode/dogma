function [f] = pair_evaluation_with_ground_truth_plot(barStruct, overlapStruct,k,iy,overlapStructGT,lengthBorders, g, ifsave)
    %   pair_evaluation_with_ground_truth_plot

    %   pair_evaluation_with_ground_truth_plot - makes a plot with the
    %   calculated and ground truth map in the same plot

    % 

    %   Returns:
    %       f,pos1, bar1, pos2, bar2,pcc,bar1Name
    %

    %   Example

    if isempty(overlapStructGT)
        numtiles = 1;
    else
        numtiles = 2;
    end

    if nargin < 8 
        ifsave = 0;
    end

    if iscell(barStruct)
        barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barStruct,'un',false);...
        cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end

%        import Core.get_full_overlap_score;
%     [fullscore,overlaplen,lenB, lenA, partialPCC,parLen,aM,bM,aP,bP] = get_full_overlap_score(oS(k,iy).pA,oS(k,iy).pB,...
%              oS(k,iy).bestBarStretch, oS(k,iy).or,barcodeGen([k iy])',oS(k,iy).h);


    % extract calculated parameters
    import Plot.ol_params_and_funs;
    [aBar,bBar, subBarA, subBarB, aFul, bFul, pos] = ol_params_and_funs(barStruct(k), barStruct(iy), overlapStruct(k,iy));
  
    if nargin< 7
        f=figure;
        g=tiledlayout(numtiles,2,'TileSpacing','compact');
    end
    ax1=nexttile(g,[1,2]);hold on
    title('A) de-novo pairwise overlap','Interpreter','latex')
    import Plot.ol_plot;
    ol_plot(aBar,bBar, subBarA, subBarB, aFul, bFul, pos,overlapStruct(k,iy),k,iy,-pos.B(1),g);

    xlim([min(pos.A(1),pos.B(1))-pos.B(1)  max(pos.A(2),pos.B(2))-pos.B(1)])

    set(gca,'ytick',[])
    xlabel('Position (kb)')
    
    xticksOriginal = get(gca, 'XTick');
    xtickLabelsOriginal = get(gca, 'XTickLabel');
    
    scaleFactor = 0.5; % should give more accurate based on bp/px
    % Scale the tick locations
    xticksScaled = xticksOriginal * scaleFactor;
%     xticksScaled = xticksScaled(xticksScaled< scaleFactor*(max(pos.A(2),pos.B(2))));
    
    
    % Optionally, update themax(pos.A(2),pos.B(2)) tick labels (if you want them to reflect the scaled values)
    xtickLabelsScaled = cellstr(num2str(xticksScaled'));
%     set(gca, 'XTick', xticksOriginal(1:length(xticksScaled)));

    set(gca, 'XTickLabel', xtickLabelsScaled);

  
    if ~isempty(overlapStructGT)
        import Plot.ol_params_and_funs;
        [aBar,bBar, subBarA, subBarB, aFul, bFul, pos2] = ol_params_and_funs(barStruct(k), barStruct(iy), overlapStructGT(k,iy));
       
        possDif = -pos2.B(1);
        ax2=nexttile([1,2]);hold on
        title('b) Reference-based pairwise overlap')
    
        import Plot.ol_plot;
        g = ol_plot(aBar,bBar, subBarA, subBarB, aFul, bFul, pos2,overlapStructGT(k,iy),k,iy,possDif,g);
        
        linkaxes([ax1 ax2],'xy')
    end

if ifsave
    [~,~] = mkdir('figs');
    saveas(f,'figs/fig3.png')
end


end

