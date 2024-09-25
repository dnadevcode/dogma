function [resultStruct,fig1] = get_display_results(barcodeGen, consensusStruct,comparisonStruct,theoryStruct, sets)
    % get_display_results 
    % Display results from the comparison of experiments vs theory
    %     Args:
    %         barcodeGen: Input 
    %         consensusStruct: 
    %         comparisonStruct: 
    %         theoryStruct:
    %         sets:
    
    resultStruct = []; % do we want to return this?
    
    % if it was chosen to display results
    if sets.displayResults ==0
        fig1 = figure('Visible', 'off');
    else
        fig1 = figure;
    end
    tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

    % choose markers for everything
    markers = ['o';'s';'x';'+';'d';'v'];

    % compute cummulative sum of lengths of barcodes
    if iscell(theoryStruct)
        lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));
    else
        lengthBorders = cumsum(arrayfun(@(x) theoryStruct(x).length,1:length(theoryStruct)));
    end
    % how many barcodes were compared to
    numBar = length(comparisonStruct);

    % If consensus was generated, then actual number of barcodes is one less
    if sets.genConsensus == 1
        numBar = numBar-length(consensusStruct);
    end

    % plot max corr coefs
%     subplot(2,2,1);hold on;
    nexttile;hold on;

    import CBT.Hca.UI.Helper.plot_max_coef;
    [fig1,maxcoef] = plot_max_coef(fig1,comparisonStruct, numBar, sets, markers);

    % plot best positions
%     subplot(2,2,2);hold on;
    nexttile;hold on;

    import CBT.Hca.UI.Helper.plot_best_pos;
    plot_best_pos(fig1,comparisonStruct, numBar, sets, markers,lengthBorders);

%     ax=subplot(2,2,3), hold on
    ax=nexttile, hold on
%     maxcoef(:,1) 
    if isequal(sets.comparisonMethod,'mp') || isequal(sets.comparisonMethod,'mpnan') || isequal(sets.comparisonMethod,'mpAll') || isequal(sets.comparisonMethod,'hmm')
        import CBT.Hca.UI.Helper.plot_best_bar_mp;
        resultStruct=plot_best_bar_mp(ax,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,1,sets);
%         [resultStruct] = plot_bar(ax, comparisonStruct, theoryStruct{1}.filename, barcodeGen, 1, 1, w, sets)
    else
        %todo: improve this plot with more information
        import CBT.Hca.UI.Helper.plot_best_bar;
        plot_best_bar(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,sets.userDefinedSeqCushion);
    end
%     end

    

            
%     subplot(2,2,4), hold on
% %     option: alternatively plot concentirc plot of the two here, based on
%     % user input
%     import CBT.Hca.UI.Helper.plot_best_image;
%     plot_best_image(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef);
% %     fig1=figure;
%         hAxis = subplot(2,2,4); hold on
    hAxis = nexttile; hold on

	if isequal(sets.comparisonMethod,'mp') || isequal(sets.comparisonMethod,'mpnan')|| isequal(sets.comparisonMethod,'mpAll') || isequal(sets.comparisonMethod,'hmm')
        sets.A = 'b';
        sets.B = 'b';
        sets.theory.isLinearTF = 1;
        resultStruct.bar1 = (resultStruct.bar1-nanmean(resultStruct.bar1))/nanstd(resultStruct.bar1);
        resultStruct.bar2 = (resultStruct.bar2-nanmean(resultStruct.bar2))/nanstd(resultStruct.bar2);
        import CBT.Hca.UI.Helper.plot_concetric;
        plot_concetric(hAxis,resultStruct,sets);
    else
        try
            import CBT.Hca.UI.Helper.plot_best_concentric_image;
            plot_best_concentric_image(hAxis,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,sets);
        catch
            % do nothing
        end
    end
    
    % option: plot all possible matches, i.e. for all 
    disp( strcat(['Number of timeframes for the unfiltered barcodes were = ' num2str(sets.timeFramesNr)]));
    
    
    
    %% additional things/options, added in HCA 4.1

    % take a simple example for concentric plot
%     idx = 87;
  
        
    % todo: make more user friendly..
    try
        mkdir(sets.output.matDirpath,sets.timestamp);
    end
    
    try
        saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'result_plot.eps'),'epsc');
    end
       
    try
        if sets.plotallmatches == 1
             mkdir(fullfile(sets.output.matDirpath,sets.timestamp),'Plots');
            for i=1:size(maxcoef,1)
                max2 = nan(size(maxcoef));      max2(i,1) = maxcoef(i,1);
                fig1 = figure('Visible', 'off');
%                 fig1 = figure;
                if isequal(sets.comparisonMethod,'mp') || isequal(sets.comparisonMethod,'mpnan')  || isequal(sets.comparisonMethod,'mpAll') || isequal(sets.comparisonMethod,'hmm')
                    ax1 = subplot(1,1,1);
                    if max2~=0
                        plot_best_bar_mp(ax1,barcodeGen,[],comparisonStruct, theoryStruct, max2,1,sets);
                    end
                else
                    plot_best_bar(fig1,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, max2);
                end
                % mp_based_on_output_pcc_test
                saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'Plots',strcat([num2str(i) '_plot.jpg'])));

%                 saveas(fig1,fullfile(sets.output.matDirpath,sets.timestamp,'Plots',strcat([sets.timestamp '_' num2str(i) '_plot.eps'])),'epsc');

            end
        end
    end
    %    assignin('base','hcaSessionStruct',hcaSessionStruct)
    
  %  cache('hcaSessionStruct') = hcaSessionStruct ;
%     end
end
