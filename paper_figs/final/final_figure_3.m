function [] = final_figure_3(sets, resRun, bG,snvr)
% Figure 3
% snr = 0.2:0.2:2;
if nargin < 4
    snrv = 1:1:10;
end
%%
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');


sets.scDiffSetting = 0.05; % scaling factor difference
sets.pxDifSetting = 50; % position difference
% sets.minOverlap = minOverlap;



% alphanu = 0.14;
alphaNu1 = 0.09;
alphaNu2 = 0.085;
alphaN1 = 0.42;
alphaN2 = 0.004;

pval = 0.05;



nBvec=[ 100 150 175 200];
setsToRun = [2 4 6 7 8 10];
% nB = 50;
ii=1;

G =[];
for ii=1:length(nBvec);
    for idxRun = setsToRun;
        nB = nBvec(ii);

        oS = resRun{idxRun}.oS;
        barcodeGen = bG{idxRun}';

        oS = oS(1:nB,1:nB);
        barcodeGen = barcodeGen(1:nB);

        % both on local and global
        [sortedVals, sortedIds,localScore,partialScore,lenA,lenB,partialLength,pvalLocal,pvalLeftOver] = ...
            calculate_sorted_pvals(oS,sets.minOverlap, alphaNu1,pval,[],alphaN1,alphaNu2,alphaN2); %localCCStruct


        sortedValsGood = sortedVals(~isnan(sortedVals));
        sortedIdsGood  = sortedIds(~isnan(sortedVals));
        localScore = localScore(~isnan(sortedVals));



        import Core.barcode_island_output;
        [barsetGen, outConsensus, coverage, consensus, islandElts, islandPx,cGenAll,barcodeIslandsData, barStruct,barIslands] =...
            barcode_island_output(sortedValsGood,sortedIdsGood, oS, barcodeGen,timestamp,sets.scDiffSetting,sets.pxDifSetting, [],0)
        %     nonEmpty= find(cellfun(@(x) length(x.members)>2,    barsetGen.barcodeIslands)); % keep only non-empty

        % nonEmpty= find(cellfun(@(x) ~isempty(x),    barsetGen.consistencyCheck)); % keep only non-empty

        import Core.create_graph_from_barsetgen;
        [Gtemp,GtempOrig] =   create_graph_from_barsetgen(barsetGen,barcodeGen,[],0);

        totalDegree = indegree(GtempOrig) + outdegree(GtempOrig);
        % Find nodes with total degree less than 2
        nodesToRemove = find(totalDegree < 2);
        G_removed = rmnode(GtempOrig, nodesToRemove);

        allNodes = 1:numnodes(GtempOrig);
        allNodes(nodesToRemove) = [];
        %     indicesCell = arrayfun(@(x) find(allNodes == x),  barIslands{idx}, 'UniformOutput', true);

        G{ii,idxRun} = G_removed;
    end
end

%%
% f = figure;
f=figure;%('Position', [10 10 500 500]);

t = tiledlayout(4,6,'TileSpacing','compact','Padding','compact')

for ii=1:length(nBvec)
    for idxRun = setsToRun
        nexttile

        nB = nBvec(ii);


        % figure
        h=plot(G{ii,idxRun},'Layout','force','ArrowSize',5,'MarkerSize',1,'NodeLabel', {});
        %

        %     highlight(h, indicesCell, 'NodeColor', 'r');
        if ii==1
            title(['', num2str(snvr(idxRun))],'Interpreter','latex')
        end

        ax = gca;

        % Rotate y-axis labels vertically
        set(ax, 'YTickLabelRotation', 90);
        if idxRun==setsToRun(1)
            ylabel(['',num2str(nB)],'FontSize',6);
        end


        if ii==3 && idxRun == setsToRun(4)
            ax2 = gca
            ax2.XTick = [];
            ax2.YTick = [];
            
            % Set other properties on the axes
            ax2.Box = 'on';
            ax2.XAxis.Color = [0 0.4470 0.7410];
            ax2.YAxis.Color = [0 0.4470 0.7410];
            ax2.LineWidth = 4;

        end

        if ii==4 && idxRun == setsToRun(5)
  
            ax2 = gca
            ax2.XTick = [];
            ax2.YTick = [];
            ax2.LineWidth = 4;
            
            % Set other properties on the axes
            ax2.Box = 'on';
            ax2.XAxis.Color = [0.8500 0.3250 0.0980];
            ax2.YAxis.Color = [0.8500 0.3250 0.0980];
  
        end

    end
end

xlabel(t,'Synthetic-noise variance ratio')
ylabel(t,'#Barcodes')



print('FIGS/Fig3.eps','-depsc','-r500');

%%

end

