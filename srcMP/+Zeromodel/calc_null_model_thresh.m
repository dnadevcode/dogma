function [mpMaxLenBasedAll] = calc_null_model_thresh(nmbpvals,fastas,theoryIdxs,dataSetIdx,nullModel,nmPerPx,numWorkers,minLen,theoryStructAll,savedata,plotdata)
    % calculating the null model:
    %
    % For hierarchical approach (Erik T. code), the estimated parameters
    % are   
    % nu = 0.198;
    % lambda = 2 * (shortLength + longLength - 2 * minOverlap); % round(0.46 * 3088286401 / 379); %
    %
    %
    import Zeromodel.ccthresh_for_comparisonstruct;

    mpMaxLenBasedAll =  cell(1,length(nmbpvals));
    switch nullModel
        case 'mpmax' % null model that calculates data self-similarity
            for nmIdx =1:length(nmbpvals)
                [~,mpMaxLenBasedAll{nmIdx},~,~,~] = bargrouping_minimum_length(fastas(theoryIdxs{dataSetIdx}),nmPerPx,nmbpvals(nmIdx),1,numWorkers,minLen, []);
            end
        case 'pval' % we calculate
            nuF = 0.12; % for nmbpvals = 0.3 check pval_test calculation for more accurate value
%             lenT = 10000;
            for nmIdx =1:length(nmbpvals)
                mpMaxLenBasedAll{nmIdx} = ccthresh_for_comparisonstruct(minLen,nuF,nuF*theoryStructAll{nmIdx}.length, 10^-4);
            end
    
        otherwise
    end
    
    if ~isempty(savedata)&& savedata==1
        save([num2str(dataSetIdx),'_','mpMaxLenBasedAll.mat'] ,'mpMaxLenBasedAll','minLen','theoryStructAll' )
    end
    
    if ~isempty(plotdata)&& plotdata==1
        
        f = figure;
        tiledlayout(2,2,'TileSpacing','compact')
        nexttile
        h = heatmap(minLen,nmbpvals,cell2mat(mpMaxLenBasedAll'));
        xlabel('Barcode length (px)');
        ylabel('nm/bp extension factor');
        h.NodeChildren(3).Title.Interpreter = 'latex';
        h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        h.NodeChildren(3).XAxis.TickLabelRotation = 60;
        % h.NodeChildren(3).XAxis.Categories(setdiff(1:length(h.NodeChildren(3).XAxis.TickValues),1:6:length(h.NodeChildren(3).XAxis.TickValues))) = [];
        title('Heatmap for $C_{thresh}$')
        
        ax = gca;
        tmp = ax.XDisplayLabels;
        tmp(setdiff(1:end,1:6:end)) = {''};
        ax.XDisplayLabels = tmp;
        
        nexttile
        nmIdx = 5;
        plot(minLen,mpMaxLenBasedAll{nmIdx})
        ylabel('$C_{thresh}$','Interpreter','latex');
        xlabel('Barcode length (px)');
        
        nexttile
        % plot comparison of distributions in case of pval:
        if isequal(nullModel,'pval')
            % compare bG to barLong
            import Nullmodel.gen_scores_test;
            [pccs,lengths,parameters] = gen_scores_test() 
%             figure
            histogram(pccs,'Normalization','pdf');
            hold on

            evp(1) = lengths(3)*0.12;
%             evp(1) = lengths(1)*0.12;
            evp(2) = 2*max(lengths)*4; % ?


            pdfF = @(cc,evdPar) evdPar(2)*(1/2*(1+betainc(cc.^2,1/2, evdPar(1)/2-1))).^(evdPar(2)-1).*(1-cc.^2).^((evdPar(1)-4)/2)./beta(1/2,  evdPar(1)/2-1);

                
            [histAll,vals] = histcounts(pccs,'Normalization','count');
            bincenters = (vals(2:end)+vals(1:end-1))/2;
            cc = bincenters;
            pdfPlot = arrayfun(@(x) pdfF(x, evp),cc);
            
            plot(cc,pdfPlot)

            xlabel('PCC score');
            ylabel('PDF')
            lgd=legend({'Histogram','Estimated PDF'})
            lgd.Location = 'southoutside';

        end

%         nexttile

    
        print('FIGS/FigS1.eps','-depsc','-r300');
    end
    
    % xlabel('Overlap length','Interpreter','latex');
    % ylabel('nm/bp extension factor','Interpreter','latex')

end

