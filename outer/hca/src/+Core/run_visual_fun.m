function [tsAlignmentVisual] = run_visual_fun(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, hcaSets,tsAlignmentVisual)

%     tsAlignmentVisual = uitab(tsAlignment, 'title', 'Visual results');


    h = tiledlayout(tsAlignmentVisual, 2,2,'TileSpacing','tight','Padding','tight');

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
    if hcaSets.genConsensus == 1
        numBar = numBar-length(consensusStruct);
    end

    % plot max corr coefs
%     subplot(2,2,1);hold on;
     ax = nexttile(h);hold(ax, 'on');grid(ax, 'on');

    import CBT.Hca.UI.Helper.plot_max_coef;
    [ax,maxcoef] = plot_max_coef(ax,comparisonStruct, numBar, hcaSets, markers);


    % plot best positions
%     subplot(2,2,2);hold on;
    ax2 = nexttile(h);hold(ax2, 'on');grid(ax2, 'on');

    import CBT.Hca.UI.Helper.plot_best_pos;
    plot_best_pos(ax2,comparisonStruct, numBar, hcaSets, markers,lengthBorders);

%     ax=subplot(2,2,3), hold on
    ax=nexttile(h); hold(ax, 'on');grid(ax, 'on');
%     maxcoef(:,1) 
    if isequal(hcaSets.comparisonMethod,'mp') || isequal(hcaSets.comparisonMethod,'mpnan') || isequal(hcaSets.comparisonMethod,'mpAll') || isequal(hcaSets.comparisonMethod,'hmm')
        import CBT.Hca.UI.Helper.plot_best_bar_mp;
        resultStruct=plot_best_bar_mp(ax,barcodeGen,consensusStruct,comparisonStruct, theoryStruct, maxcoef,1,hcaSets);
%         [resultStruct] = plot_bar(ax, comparisonStruct, theoryStruct{1}.filename, barcodeGen, 1, 1, w, hcaSets)
    else
        %todo: improve this plot with more information
        import CBT.Hca.UI.Helper.plot_best_bar;
        plot_best_bar(ax,barcodeGenC,consensusStruct,comparisonStruct, theoryStruct, maxcoef,hcaSets.userDefinedSeqCushion);
    end

    hAxis = nexttile(h); hold(ax, 'on');grid(ax, 'on');
    import CBT.Hca.UI.Helper.plot_length_vs_intensity;
    plot_length_vs_intensity(hAxis,barcodeGenC)

    % option: plot all possible matches, i.e. for all 
    disp( strcat(['Number of timeframes for the unfiltered barcodes were = ' num2str(hcaSets.timeFramesNr)]));
    
    %%
    if hcaSets.plotallmatches == 1 
        disp( strcat(['Started saving plots ..']));
        tic;
        if hcaSets.plotalldiscrim && isfield(comparisonStruct{1},'discriminative')
            idxToPlot = find(cellfun(@(x) x.discriminative.is_distinct,comparisonStruct));
        else
            idxToPlot = 1:size(maxcoef,1);
        end

        [~,~] = mkdir(hcaSets.output.matDirpath,'all_molecules');
        [~,~]= mkdir(fullfile(hcaSets.output.matDirpath,'all_molecules',[hcaSets.timestamp,'Plots']));

        saveas(h,fullfile(hcaSets.output.matDirpath,'all_molecules',[hcaSets.timestamp,'Plots'],'result_plot.fig'));
    
        for i=idxToPlot
%             tic
                max2 = nan(size(maxcoef));      max2(i,1) = maxcoef(i,1);
                fig1 = figure('Visible', 'off');
                ax1 = subplot(fig1,1,1,1);


%                 fig1 = figure('Visible', 'on');

                if isequal(hcaSets.comparisonMethod,'mp') || isequal(hcaSets.comparisonMethod,'mpnan')  || isequal(hcaSets.comparisonMethod,'mpAll') || isequal(hcaSets.comparisonMethod,'hmm')
                    if max2~=0
                        plot_best_bar_mp(ax1,barcodeGenC,[],comparisonStruct, theoryStruct, max2,1,hcaSets);
                    end
                else
                    hold on
                    plot_best_bar(ax1,barcodeGenC,consensusStruct,comparisonStruct, theoryStruct, max2);
                end

                saveas(fig1,fullfile(hcaSets.output.matDirpath,'all_molecules',[hcaSets.timestamp,'Plots'],strcat([num2str(i) '_plot.jpg'])));
                close(fig1);
%                 saveas(fig1,fullfile(hcaSets.output.matDirpath,hcaSets.timestamp,'Plots',strcat([hcaSets.timestamp '_' num2str(i) '_plot.eps'])),'epsc');
        
        end
        timePassed = toc;
        disp(strcat(['Experiments vs theory plots were saved in ' num2str(timePassed) ' seconds']));

    end
    %%

end

