% Evaluation figure - experiments from individual days
%
% Individual barcodes can be compare against the theory, using
%
% 1) Standard PCC (+ pvalue to sort)
% 2) MP (+ p-value to sort)
% 3) Full discriminative barcode analysis against whole database
%
% As output, gives a position and significance

%% PARAMETERS
% psffac = 1; % scaling factor (in case PSF needs to be something else than 300nm)
numWorkers = 30; % num workers, in one node there are 30
minLen = [150:50:3000]; % min to max % could take more points for more accurate..
sets.comparisonMethod = 'mass_pcc';
sF = 0.9:0.01:1.1;
nmPerPx = 110;

% choose null-model
nullModel = 'mpmax'; % 'pval'
nullModel='pval';

%% Gen theories all possible nmbp range
% if for all nm/bp values
nmbpvals = 0.15:0.01:0.3; % should use mpMax to get this "fully" correct
import Thry.gen_theoretical;
theoryStructAll = cell(1,length(nmbpvals));
for nmIdx =1:length(nmbpvals)
    [theoryStructAll{nmIdx},~,barcodeGenT] = gen_theoretical(fastas(theoryIdxs{dataSetIdx}),nmbpvals(nmIdx),0,kymoStructs{1}{1}.nmpxidFold); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly
end

%% Null model, calculates threshold for PCC score
import Zeromodel.calc_null_model_thresh;
[mpMaxLenBasedAll] = calc_null_model_thresh(nmbpvals,fastas,...
    theoryIdxs,dataSetIdx,nullModel,nmPerPx,numWorkers,minLen,theoryStructAll,1,0);

mpMaxLenBasedC = cell(1,length(bG));
mpMaxLenBasedClocal = cell(1,length(bG));
thrS = cell(1,length(bG));


for idxRun = 1:length(bG)
    idxRun
    try
    nmPerPx = kymoStructs{idxRun}{1}.nmpxidFold;
    nmbp = kymoStructs{idxRun}{1}.nmBpidFold;
    catch
        warning('params not provided')
        nmPerPx = 110;
        nmbp = 0.25;
    end
    fastaFile = fastas((theoryIdxs{dataSetIdx}));
%     [mpMaxLenBasedClocal{idxRun}] = calc_null_model_thresh(nmbp,fastas,...
%     theoryIdxs,dataSetIdx,'mpmax',nmPerPx,numWorkers,minLen,theoryStructAll,1,0);

  %     figure,plot(mpMaxLenBasedClocal{idxRun}{1});hold on;plot(mpMaxLenBasedC{idxRun}{1})

    [thrS{idxRun},~,barcodeGenT] = gen_theoretical(fastas(theoryIdxs{dataSetIdx}),nmbp,0,kymoStructs{1}{1}.nmpxidFold); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly
%     [MP,mpMaxLenBasedC{idxRun},theoryStructRev,MPI,thrS{idxRun}] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,1,numWorkers,minLen, []);% no lenseq if based on fastafile
  [mpMaxLenBasedC{idxRun}] = calc_null_model_thresh(nmbp,fastas,...
    theoryIdxs,dataSetIdx,nullModel,nmPerPx,numWorkers,minLen,{thrS{idxRun}},1,0);

end

calcNMrate = 1;
plotExample = 0;

successRate = zeros(1,length(bG));
sucRateStruct = cell(1,length(bG));
comparisonStructC= cell(1,length(bG));
passthreshC= cell(1,length(bG));
rezMaxAll = cell(1,length(bG));

for idxRun = 1:length(bG)
    idxRun

    % compare to theory: either pcc or mp
    [comparisonStructC{idxRun},rezMax,bestBarStretch] = compare_to_t(bG{idxRun},thrS{idxRun},sF,sets);

    % comparison struct to overlap struct:
    import Core.comparisonsres_to_struct; 
    [oSreal] = comparisonsres_to_struct(bG{idxRun}, comparisonStructC{idxRun},thrS{idxRun}.length);


    allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStructC{idxRun});
    allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStructC{idxRun});
    
%         [mpMaxLenBasedC{idxRun}] = calc_null_model_thresh(nmbp,fastas,...
%     theoryIdxs,dataSetIdx,nullModel,nmPerPx,numWorkers,minLen,thrS{idxRun},1,0);

    idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
    coefDiffs = allCoefs-mpMaxLenBasedC{idxRun}{1}(idxThresh);
    
    passthreshC{idxRun} = find(coefDiffs > 0);
    successRate(idxRun) = length(passthreshC{idxRun} )/length(coefDiffs); % this could still have false positives.. ?
    
    % alternative
    pasthresh = find(cellfun(@(x) x.pval, comparisonStructC{idxRun}) < 0.0001);

    if plotExample

%         pvals = pvalfun(allCoefs,allLengths)
        % fig for paper?
        import CBT.Hca.UI.Helper.plot_any_bar;
        plot_any_bar(1,bG{idxRun},{comparisonStructC{idxRun}},thrS{idxRun},1);
    end

    % calc p-values /full

    if calcNMrate
% %     % check for different nmbp (in case estimated incorrectly
    nmbpvals = 0.15:0.01:0.3; % should use mpMax to get this "fully" correct
    import Thry.gen_theoretical;
    succesNmnp = zeros(1,length(nmbpvals));
    coefDiffsAll = zeros(1,length(nmbpvals));

    for nmIdx =1:length(nmbpvals)
        %         [theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbpvals(nmIdx),0,nmPerPx); %todo faster: length rscale instead of calculating each time. also this can be done in outside loop. also make it into struct so don't need to load it here and manipulate theory directly
        [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(bG{idxRun},theoryStructAll{nmIdx},sF,sets);
        
        allCoefs = cellfun(@(x) x.maxcoef(1),comparisonStruct);
        allLengths = cellfun(@(x) x.lengthMatch (1),comparisonStruct);
        
        idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
        coefDiffs = allCoefs-mpMaxLenBasedAll{nmIdx}(idxThresh);
        
        barsPassThresh = find(coefDiffs > 0.05);
        succesNmnp(nmIdx) = length(barsPassThresh)/length(coefDiffs); % this could still have false positives.. ?
%         succesNmnp(nmIdx) = length(find(cellfun(@(x) x.maxcoef(1),comparisonStruct)-mpMaxLenBasedAll{nmIdx}(idxThresh) > 0.05))/length(bG{idxRun}); % this could still have false positives.. ?
        coefDiffsAll(nmIdx) = mean(coefDiffs);
    % todo: same with MP
    end
    sucRateStruct{idxRun}.succesNmnp = succesNmnp;
% 
    f = figure
    plot(nmbpvals,sucRateStruct{idxRun}.succesNmnp );
    xlabel('nm/bp','Interpreter','latex')
    ylabel('succes rate','Interpreter','latex');
    hold on
    plot(nmbp,0,'redx')
    lgd = legend({'Success rate','Estimated nm/bp'})
    lgd.Location ='southoutside';
    title('Success rate for different nm/bp extension factors','Interpreter','latex')
%     print('FIGS/FigS2.eps','-depsc','-r300');
% 
    end
end

save([num2str(dataSetIdx),'_individualdays.mat'] ,'comparisonStructC','thrS','bG','sets','sucRateStruct','mpMaxLenBasedC')
%%
idxRun = 1;
f = figure
plot(nmbpvals,   sucRateStruct{idxRun}.succesNmnp);
xlabel('nm/bp','Interpreter','latex')
ylabel('succes rate','Interpreter','latex');
hold on
plot(kymoStructs{idxRun}{1}.nmBpidFold,0,'redx')
lgd = legend({'Success rate','Estimated nm/bp'})
lgd.Location ='southoutside';
title('Success rate for different nm/bp extension factors','Interpreter','latex')
%
%%
nmbp = kymoStructs{idxRun}{1}.nmBpidFold;
goodbars = bG{idxRun}(1);
w = [0 400]; % only global
import Core.load_theory_structure;
thryFileIdx = 1; % todo: pass directly the theory file here
[theoryStruct,sets] = load_theory_structure(nmbp,thryFileIdx);

import Core.extract_species_name;
[speciesLevel,idc] = extract_species_name(theoryStruct);
% 
[rezOut] = local_alignment_assembly(theoryStruct,goodbars,w);

idx = 1;
import Core.disc_locs;
[refNums, allNums, bestCoefs,refNumBad, bestCoefsBad] = disc_locs(rezOut{idx}.rezMax);
{theoryStruct([refNums{1}]).name}'

% local bootstrap
[scores,pccScore] = local_bootstrap_run(goodbars,cellfun(@(x) x.rezMax,rezOut,'un',false), repmat({matlab.lang.makeValidName(strrep( goodbars{idx}.name,'.tif',''))},1,length(w)), theoryStruct ,w, speciesLevel, idc);

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
import Core.export_coefs_resampling;
T = export_coefs_resampling(scores,goodbars, w, [pwd, '/resampling_table'],timestamp);

super_quick_plot(1,goodbars,rezOut{1},theoryStruct)
super_quick_plot_rezmax(1,goodbars,rezOut{1}.rezMax,theoryStruct,refNums{1}(1))
% 
% import Core.discrim_true_positives;
% [truePositives,discSpecies,discAll,allSpecies,refNums,signMatch] =...
% discrim_true_positives(rezOut{idx}.rezMax,speciesLevel,idc);




% %%
% bpPx = 500;
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], comparisonStructC, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6
% % 
idxRun = 9
bgen =  bG{idxRun};
[outConsensus, coverage, pval] = gen_reference_based_assembly(bgen,comparisonStructC{idxRun},thrS(idxRun),'test11',inf);
% 
% 
% % quick_visual_plot(46,1,bgAll,rezMax,bestBarStretch,theoryStruct)
% 
% [compStr,~,calcLengths] = compare_to_t_mp(bgAll,theoryStruct,sF,300); %todo: extend local to full
% [comparisonStruct{4}.pos(1) compStr{2}.pos(1)]
% 
% passthreshLen = find(ismember(calcLengths,barsPassThresh));
% % passthreshLen = barsPassThresh(calcLengths);
% 
% bpPx = 500;
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], compStr(passthreshLen), [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6
% 
% % compStr{passthreshLen(1)}.pos = 4169;
% % super_quick_plot(1,bgAll(4),compStr(2),theoryStruct)
% 
% % allCoefs = cellfun(@(x) x.maxcoef(1),compStr);
% % allLengths =300;
% % 
% % idxThresh = arrayfun(@(x) find(x<minLen,1,'first'),allLengths);
% % coefDiffs = allCoefs-mpMaxLenBased(idxThresh);
% % 
% % barsPassThresh = find(coefDiffs > 0);
% % 
% % bg = bgAll(barsPassThresh); % make gen_reference_based_assembly it work for MP
% [outConsensus, coverage, pval] = gen_reference_based_assembly(bgAll(passthreshLen),compStr(passthreshLen),theoryStruct,'test11',inf);
% 
