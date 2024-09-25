function [outC] = plot_island(posStruct, islandElts, barcodeGen, toplot,tsVisual,    numStds)
%   Creates islands

    if nargin < 5 || isempty(tsVisual)
    
        %   TODO: plot in the tabbedscreen
          % create tabbed figure
        hFig = figure('Name', 'Barcode Islands', ...
            'Units', 'normalized', ...
            'OuterPosition', [0 0 0.9,0.9], ...
            'NumberTitle', 'off', ...
            'MenuBar', 'none', ...
            'ToolBar', 'figure' ...
        );
        hPanel = uipanel('Parent', hFig);
        h = uitabgroup('Parent',hPanel);
        t1 = uitab(h, 'title', 'islands');
       
        tsHCC = uitabgroup('Parent',t1);
    else
        tsHCC = uitabgroup('Parent',tsVisual,'Units', 'normalized','Position', [0.01 0.01 0.99 0.99]);
    end

    if nargin < 6
        numStds = 2;
    end

    ordIdx = cellfun(@(x) x.barid,posStruct); % order of barcodes

    % convert the reference coefs to p-vals, etc.
    orientation = cell2mat(cellfun(@(x) x.or,posStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,posStruct,'UniformOutput',0)');
    maxcoef = cell2mat(cellfun(@(x) x.maxcoef,posStruct,'UniformOutput',0)');
    strF = cell2mat(cellfun(@(x) x.bestBarStretch,posStruct,'UniformOutput',0)');

    lengthsbar = cellfun(@(x) length(x.rawBarcode),barcodeGen);

    outC = cell(1,length(islandElts));
    t =[];
    for  i=1:length(islandElts)
        posI = pos(islandElts{i});
        orientationI = orientation(islandElts{i});
        maxcoefI = maxcoef(islandElts{i});
        strFI = strF(islandElts{i});
        lengthsbarI = lengthsbar(islandElts{i});


        posI = posI-min(posI)+1;   
        
        

        %% plot multi theories with some nan's in between
        maxLen = max(posI+round(lengthsbarI.*strFI));
        outConsensus = nan(length(lengthsbarI),maxLen);
    %     outConsensus(1,:) = zscore(theorBar);
%     lenT = length(theorBar);
    
    for k=1:length(lengthsbarI)
%         ii = goodBars(k);
        expBar = barcodeGen{islandElts{i}(k)}.rawBarcode;
        expBit = barcodeGen{islandElts{i}(k)}.rawBitmask;

%         expLen = length(expBar);

        % interpolate to the length which gave best CC value
        expBar = imresize(expBar,'Scale' ,[1 strFI(k)]);
        expBit = imresize(expBit,'Scale' ,[1 strFI(k)]);

        expBar(~expBit)= nan;

        if orientationI(k,1)==2
            expBar = fliplr(expBar);
            expBit = fliplr(expBit);
        end
        outConsensus(k,posI(k):posI(k)+length(expBar)-1) = (expBar-nanmean(expBar))/nanstd(expBar,1)+k*0;

    end

    % test: remove pixels which are outliers
    lowerBound = nanmean(outConsensus)-numStds*nanmean(nanstd(outConsensus));
    higherBound = nanmean(outConsensus)+numStds*nanmean(nanstd(outConsensus));

    outConsensus(outConsensus<lowerBound) = nan;
    outConsensus(outConsensus>higherBound) = nan;


    outC{i} = outConsensus;

% h1 = uitabgroup('Parent',hPanelPlot);
%    hTabgroup = uitabgroup('Parent',hPanelResult{i});
        hPanelPlot = uitab(tsHCC, 'title', strcat(['Bar Island_' num2str(i)]));

%         hResScores= uitab(hTabgroup, 'title',strcat('Scores'));
        t{i} = tiledlayout(hPanelPlot,2,1,'TileSpacing','tight','Padding','tight');
      
%     figure,
%         tiledlayout(2,1);
        nexttile(t{i})
        plot(outConsensus')
        nexttile(t{i})
        plot(outConsensus'+[1:size(outConsensus,1)]*3)
        legend({num2str((ordIdx(islandElts{i})))})
        %arrayfun(@(x) num2str(x), (ordIdx(islandElts{i})), 'UniformOutput', false)
    
    
    end
    
%     ccthresh = 0.3;
%     nuF = 0.03;
%     import Zeromodel.beta_ev_cdf; % do we need here? already best barcodes by sorted values..
%     pval = nan(1,length(barcodeGen));
%     for i=1:length(barcodeGen)
% 
%         if comparisonStruct{i}.maxcoef > ccthresh % this is for overlap, so change..
%             pval(i) = 1-beta_ev_cdf(comparisonStruct{i}.maxcoef, nuF*lengthsbar(i), 1, 2*max(lengthsbar(i),lengthsThry(idxthry(i))),1);
%         end
%     end

% %     for i=1:length(barcodeGen)
% %         min_overlap = lengthsbar(i); %min(lengthsThry(i),lengthsbar(i));
% %         a = 0.13*min_overlap;
% %         n = 2*(lengthsThry(i)-min_overlap);
% %         pval(i) = 1-Zeromodel.beta_ev_cdf(maxcoef(i), a, 1, n, 1);
% % 
% % %         0.004*(len2-MIN_OVERLAP_PIXELS+1)*(len1-MIN_OVERLAP_PIXELS+1);
% %     end
    
%     goodBars = find(pval< 0.01);
    
    % load theory file. just single one
%     fileID = fopen(theoryStruct{comparisonStruct{1}.idx}.filename,'r');
%     formatSpec = '%f';
%     theorBar = fscanf(fileID,formatSpec);
%     fclose(fileID);
  
%     
% PCC_OVERLAP = reshape([overlapStruct.fullscore], size(overlapStruct,1),size(overlapStruct,2));
% lenRoot = reshape([overlapStruct.lenB], size(overlapStruct,1),size(overlapStruct,2));
% lenA = reshape([overlapStruct.lenA], size(overlapStruct,1),size(overlapStruct,2));
% overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));
% 
% % index of ROOT barcode 
% pccVals = PCC_OVERLAP(:,root); % columns are non re-scaled
% 
% % overlapStruct
% %%
% % bar1 = 20;
% % figure,plot(PCC_MP(bar1,:))
% % figure,plot(PCC_OVERLAP(bar1,:))
% 
% idxToCheck = find(pccVals>0.5); % check only id's passing initial thresh
% 
% pval = nan(1,length(pccVals));
% import Zeromodel.beta_ev_pdf;
% for j=idxToCheck'
%         a = 0.08*overlaplen(j,root); % based on pval_test tests
%         n = 2*(lenRoot(j, root)+lenA(j,root)-2*overlaplen(j,root));
%         pval(j) = 1-Zeromodel.beta_ev_cdf( PCC_OVERLAP(j,root), a, 1, n, 1);
% end
% 
% numGoodMatch = pval < pthresh;
% 
% goodvals = pval(pval < pthresh);
% peaksToTry = find(numGoodMatch);
% 
% orBars = [overlapStruct(peaksToTry,root).or];
% reFac = [overlapStruct(peaksToTry,root).bestBarStretch];
% pA = [overlapStruct(peaksToTry,root).pA];
% pB = [overlapStruct(peaksToTry,root).pB];
% barsFinal=[root peaksToTry];
% 
% if ~isempty(toplot)
%     fig = figure; 
%     tiledlayout(1,2);
%     nexttile([1 2])
%     hold on
%     plot(zscore(barStruct(root).rawBarcode(barStruct(root).rawBitmask)),'black')
%     for j=1:length(peaksToTry)
%         rescaled = imresize(barStruct(peaksToTry(j)).rawBarcode(barStruct(peaksToTry(j)).rawBitmask),'Scale' ,[1 reFac(j)]);
%         if orBars(j)==-1
%             rescaled = fliplr(rescaled);
%         end
%         plot(pB(j)-pA(j)+1:pB(j)-pA(j)+length(rescaled),zscore(rescaled)+3*j,'black')
%     end
% 
%     bars = cell(1,length(peaksToTry)+1);
%     bars{1} = strcat('root ',num2str(root));
% %     [num2str(root) arrayfun(@(x) num2str(x),[peaksToTry],'un',false)];
%     for j=1:length(peaksToTry)
%          bars{j+1} = strcat([num2str(peaksToTry(j)) '; or=' num2str( orBars(j))...
%              '; rf=' num2str(reFac(j)) '; pcc=' num2str(pccVals(peaksToTry(j)),'%3.2f') '; pval=' num2str(goodvals(j),'%4.2e')])
%     end
%         legend(fliplr(bars),'location','southoutside','Interpreter','latex');
%         saveas(fig,'figs/fig6.png')
% 
% %      for j=1:length(peaksToTry)+1
% %          plot( )
% 
% %     if nargin>=12
% %         saveas(f,'figs/bargroupexample.png')
% %     end

end
% 
%     % bargroup..
%     import Core.plot_bargroup;
%     if pl==1
%         [barMat,bars,orBars,reFac] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,pl);
%     else
%         [barMat,bars,orBars,reFac] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
%     end        
%     % instead create matrix like in gen reference based assembly
% %     figure;hold on;
% %     for i=1:length(barMat)
% %         plot(  barMat{i}{1},zscore(  barMat{i}{2} )+3*i,'red');
% %     end

% end

