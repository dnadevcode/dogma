function [] = dogma_run(bargiSets, pathMain, userDir,tsVisual)
    % bargi_run - BARGrouping Islands run script
    %
    %   Example:
    %       
    %   [bargiSets] = Core.Default.read_default_bargi_sets('bargi_settings.txt');
    %   bargi_run(bargiSets);
    %     bargi_run([],'/home/avesta/albertas/reps',foldtests)

   

    %   Args:
    %       pathMain - path where the required repositories are
    %       dbmOSW - settings structure

    %   Returns:

    if nargin < 1 || isempty(bargiSets)
         [bargiSets] = Core.Default.read_default_bargi_sets('bargi_settings.txt');
         bargiSets = bargiSets.default;
    end   
    if nargin < 3 && ~isfield(bargiSets,'kymofolder') && ~isequal(bargiSets.default.method,'synth')
        bargiSets.kymofolder = uigetdir(pwd); 
    else
        if ~isempty(userDir)
            bargiSets.kymofolder = userDir;
        end
    end

    if nargin < 4

        mFilePath = fileparts(fileparts(mfilename('fullpath')));
        versionBargi = importdata(fullfile(mFilePath,'VERSION'));
        bargiSets.version =  versionBargi{1};
        hFig = figure('Name', ['Bardenas ' versionBargi{1}], ...
            'Units', 'normalized', ...
            'OuterPosition', [0.05 0.1 0.8 0.8], ...
            'NumberTitle', 'off', ...     
            'MenuBar', 'none',...
            'ToolBar', 'none' ...
        );
        hPanel = uipanel('Parent', hFig);
        h = uitabgroup('Parent',hPanel);
        ts = uitab(h, 'title', 'Bargi settings');
        tsVisualV = uitab(h, 'title', 'Visual results');
        tsVs = uitabgroup('Parent',tsVisualV);
        tsVisual.single = uitab(tsVs, 'title', 'Single');
        tsVisual.islands = uitab(tsVs, 'title', 'Islands');
        tsVisual.graph = uitab(tsVs, 'title', 'Graph');
        tsVisual.block = uitab(tsVs, 'title', 'Block');

    end

    if ~isfield(bargiSets,'kymofolder')
        bargiSets.kymofolder = [];
    end

    set_up(); % set up things

    % time-stamp for the results of this run
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

    sets = bargiSets;
    sF = 1-sets.sF/100:sets.sFstep/100:1+sets.sF/100; % length re-scaling factor

    % three options: sample 2, sample 5, and ecoli 2
    import Core.load_chrom_data;
    [bgAll, bG, kymoStructs, synthStr, synthStr2, theoryStruct] = load_chrom_data(bargiSets.kymofolder,sets.method,sets.depth,sets.numframes);

    if iscell(bargiSets.kymofolder)
        outfold =fileparts(bargiSets.kymofolder{1});
    else
        outfold =fileparts(bargiSets.kymofolder);
    end
    % TODO: filter her based on some infoscore
%     try
%         snr =  cellfun(@(x)  x.snr,bgAll);
%         bgAll((snr>10)) = [];
%     catch 
% 
%     end


    
    bars = bgAll(cellfun(@(x) sum(x.rawBitmask),bgAll)*sF(1)>sets.minOverlap); % only for those larger than min overlap
    tic
    [oS] = calc_overlap_mp(bars,sF, sets.minOverlap,timestamp,sets.numworkers );
    disp(['All pairwise overlaps calculated in ', num2str(toc), ' sec']);
    
    save(fullfile(outfold,[sets.method,timestamp,'_mp_data.mat']),'sets','oS','sF','bars')

    % after this part, it's only visualization, most of heavy calculation
    % is already done 
        
    [sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
        calculate_sorted_pvals(oS,sets.minOverlap, sets.nupar,sets.pthresh);
%     [sortedVals, sortedIds,pscores,fullScores,overlaplen,partialScore,partialLength,pval,pvalLeftOver,pvalCombined] = sorted_scores(oS);
    
    numSigOverlap = sum(~isnan(sortedVals));
    disp(['Found ', num2str(numSigOverlap) ,' significant overlaps'])


% plot1: take ground trtuh comparisonStruct

%     if plotfigs
%     figure,histogram(sortedVals)

idxPair = 1;
[curSink, curSource] = ind2sub(size(oS),sortedIds(idxPair));

%         lenThry = theoryStruct{1}{1}.length;
import Plot.pair_evaluation_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)

if sets.plotfigs && numSigOverlap > 0
    lenThry = 10000;
    h = tiledlayout(tsVisual.single,2,4,'TileSpacing','tight','Padding','tight');
    if ~isempty(synthStr)
        [h] = pair_evaluation_plot(bars, oS,curSink,curSource,synthStr{1}([curSink curSource]),lenThry,0,h);
    end
    pair_evaluation_with_ground_truth_plot(bars, oS,curSink,curSource,[],lenThry,h);
    if sets.savefigs
        exportgraphics(tsVisual.single,fullfile(outfold,[sets.method,timestamp,'overlap.png']))
    end
end
%     end

barcodeGen = bars;
idxRun = 1;
% oS = resRun{idxRun}.oS;
% barcodeGen = bG{idxRun};

goodPos = sortedVals > sets.zscoreCutOff;
sortedIdsGood = sortedIds(find(goodPos));
sortedValsGood = sortedVals(find(goodPos));

% import Core.create_overlap_graph;
% [finalgraph,Ggraphs] = create_overlap_graph(size(oS),sortedIdsGood,length(sortedIdsGood),bars,zeros(1,length(sortedIdsGood)),100,oS,tsVisual.graph);


% refinement to make graph circular? do we need to remove weak links (that
% might represent false mappings?)

% plot output
import Core.barcode_island_output;

[barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGen,barcodeIslandsData, barStruct,barIslands]=...
    barcode_island_output(sortedValsGood,sortedIdsGood, oS, bars,timestamp,sets.scDiffSetting,sets.pxDifSetting, tsVisual.islands);

    if sets.savefigs
        for i=1:length(tsVisual.islands.Children)
            if ~isempty(tsVisual.islands.Children(i).SelectedTab)
                exportgraphics(tsVisual.islands.Children(i).SelectedTab,fullfile(outfold,[sets.method,timestamp, num2str(i), 'islands.png']));
            end
        end
    end


    disp(['Found ', num2str(length(outConsensus)) ,' islands'])

    import Core.create_graph_from_barsetgen;
    create_graph_from_barsetgen(barsetGen,bars,tsVisual.graph);

    consistency = barsetGen.consistencyCheck(cellfun(@(x) ~isempty(x), barsetGen.consistencyCheck));

    save(fullfile(outfold,[sets.method,timestamp,'_res_data.mat']),'sets','barsetGen','outConsensus','cGen');


%     figure,plot(cellfun(@(x) x.positionalDifference,consistency))
%     figure,plot(cellfun(@(x) x.failedfilt,consistency))
%     figure,plot(cellfun(@(x) x.sFratio,consistency))

    % block representation
    if sets.plotfigs
        import Plot.plot_bloc_rep;
        plot_bloc_rep(outConsensus,tsVisual);
        if sets.savefigs
            for i=1:length(tsVisual.block.Children)
                if ~isempty(tsVisual.block.Children(i).SelectedTab)
                exportgraphics(tsVisual.block.Children(i).SelectedTab,fullfile(outfold,[sets.method,timestamp, num2str(i), 'block.png']));
                end
            end
        end
    end


    if sets.saveOutput
            assignin('base','oS', oS)
            assignin('base','bars', bars)
            assignin('base','outConsensus', outConsensus)
            assignin('base','barsetGen', barsetGen)
            assignin('base','sets', sets)
    end



    function set_up()    
        % Load required packages:
        if isempty(pathMain) || nargin < 2
            pathMain = 'C:\Users\Lenovo\git\dnadevcode\'; % path to reps.
        end    
        addpath(genpath([pathMain,'lldev'])); % should not need lldev
        addpath(genpath([pathMain,'hca']));
        addpath(genpath([pathMain,'bargroupingprototype']));
    end


end

