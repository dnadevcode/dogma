function [tpr, fpr] = bml_validation(synthStr,bG,theoryStruct,N, pxpsf, mL, stdL, ...
    snr, stdE, isC, tL, minL )

    %% min length calculation for different overlap lengths. Gives best possible score for each length
    % could also calculate sigma (bootstrapping the best length?)
    pxpsf = 110;
    nmPerPx = 300/pxpsf;
    nmbp = 0.2; % length re-scale factor

    minLen =[minL:50:3000]; % min to max
    
    nullModel = 'pval';
    import Zeromodel.calc_null_model_thresh;
%     [mpMaxLenBased] = calc_null_model_thresh(nmbpvals,...
%     theoryIdxs,dataSetIdx,nullModel,nmPerPx,numWorkers,minLen,theoryStructAll,1,0);

    [mpMaxLenBased] = calc_null_model_thresh(nmbp,fastas,...
    theoryIdxs,dataSetIdx,nullModel,nmPerPx,numWorkers,minLen,theoryStruct{1},1,0);
% nmbpvals,fastas,theoryIdxs,dataSetIdx,nullModel,nmPerPx,numWorkers,minLen,theoryStructAll,savedata,plotdata

%     [MP,mpMaxLenBased,theoryStructRev,MPI,~] = bargrouping_minimum_length([],nmPerPx,nmbp,1,30,minLen, tL*nmPerPx/nmbp);

%%
    %
    % SNR = 1:10;
    tpr = zeros(1,length(snr));
    fpr = zeros(1,length(snr));
    
    cdiff = 0.05; % cut-off above minLen to say that it is significant. Alternative: get this by bootstrapping best length
    pvalthresh = 0.01;
    sets.comparisonMethod = 'mass_pcc';
    sF = 0.9:0.025:1.1; % same as species paper
    
    
    import Zeromodel.beta_ev_cdf;
    nuF = 0.12; % should be tuned based on synthetic data? Depends a bit on pdf
    
    allowedDist = 25;% allowed distance, this can be calculated from barLength*stepE, where currently stepE = 0.025
    
            
    % maybe show this for varying allowedDist?
    
    for i=1:length(snr)
        % compare. maybe use local_alignment_assembly function here?
        [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(bG{i},theoryStruct{i},sF,sets);
    
        % get p-value
        maxcoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
        pval1 = arrayfun(@(x,y) 1-beta_ev_cdf(x, nuF*y, 1, nuF*2*max(y, theoryStruct{i}{1}.length),0),maxcoefs,...
            cellfun(@(x) x.lengthMatch,comparisonStruct));
    
        
        allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
        allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);
        
        idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
        coefDiffs = allCoefs-mpMaxLenBased{1}(idxThresh);
        
        barsParsThresh = coefDiffs > cdiff;
    %     barsParsThresh = pval1 < pvalthresh; % p-val doesn't take into
    %     account possible sequence similarity regions within the sequence
    
        % plot to show that length thresholding works similar to match score
        % thresholding?
%         figure,plot(barsParsThreshPval==0,'x');hold on;plot(barsParsThresh==0,'o')
    
        % barsPassThresh = find();
    
        
%         super_quick_plot(6,barcodeGen,comparisonStruct,theoryStruct)
        
        %% check success rate statistics
        % N = length(barcodeGen);
        posShift = zeros(1,N);
        tp = zeros(1,N);
        fn = zeros(1,N);
        fp = zeros(1,N);
        tn = zeros(1,N);
        
        success = 0;
        for idx=1:N
            pos = comparisonStruct{idx}.pos(1);
            posGT = synthStr{i}{idx}.pos;
            %         posDifTheory = pos-posGT;
            posShift(idx) = min([dist(pos,posGT) dist(pos+tL,posGT) dist(pos-tL,posGT) ]);
            if posShift(idx) <= allowedDist
                if barsParsThresh(idx)==1
                    tp(idx) = 1;
                else
                    fn(idx) = 1;
                end
            else
                if barsParsThresh(idx)==1
                    fp(idx) = 1;
                else
                    tn(idx) = 1;
                end
            end
        end
        
        tpr(i) = sum(tp)/length(bG{i});%(sum(tp)+sum(fn)); % sensitivity
        fpr(i) = sum(fp)/(sum(fp)+sum(tn)); % fprc
    end
    save(['synth_',timestamp,'bml.mat'] ,'tpr','fpr', 'comparisonStruct','theoryStruct')

    f = figure;
    tiledlayout(1,2,'TileSpacing','tight');
    % plot 1 success rate
    nexttile
    plot(tpr)
    hold on
    plot(fpr)
    xlabel('SNR')
    ylabel('Success rate')
    legend({'tpr','fpr'},'location','southoutside')
    title('(A) Validation of ground truth method','Interpreter','latex')
    % plot 2 match positions against theory
    nexttile
%     title
    import Plot.plot_map_on_reference;
    f = plot_map_on_reference(f, comparisonStruct, [], theoryStruct{i}{1}.length,1);%1/bpPx*10^6
    title('(B) Positions along the reference','Interpreter','latex')

    print('FIGS/FigS10.eps','-depsc','-r300');


%     [outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,comparisonStruct{idxRun},thrS(idxRun),'test11',inf);




end

