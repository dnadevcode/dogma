function [comparisonStruct,rezMax,bestBarStretch] = compare_to_t(barcodeGen,theoryStruct,sF,sets)
    % to remove/same as old fun in hca
%% PCC
% barcodeGen = barcodeGen1;
if isempty(sF)
    sF = 0.8:0.01:1.2;
end

sets.filterSettings.filter = 0;

if ~isfield(sets,'w')
    sets.w = [];
end
rezMax=[];bestBarStretch=[];bestLength=[];
for i=1:length(theoryStruct)
%     tic
    import CBT.Hca.Core.Comparison.on_compare;
    if iscell(theoryStruct)
    [rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGen,...
        theoryStruct{i},sets.comparisonMethod, sF,sets.w,50,[],sets.filterSettings);
    else
        [rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGen,...
        theoryStruct(i),sets.comparisonMethod, sF,sets.w,50,[],sets.filterSettings);

    end
%     toc
end

comparisonStructAll = rezMax;
for i=1:length(comparisonStructAll)
    for j=1:length(bestBarStretch{i})
        comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
        comparisonStructAll{i}{j}.length = bestLength{i}(j);
    end
end
import CBT.Hca.Core.Comparison.combine_theory_results;
[comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);

% pval:
if ~isfield(sets,'nuF')
    alphaNu = 0.12; % depends on PSF, here fix for 110nm/px
else
    alphaNu = sets.nuF;
end

if ~isfield(sets,'nF')
    alphaN = 0.12; % depends on PSF, here fix for 110nm/px
else
    alphaN = sets.nF;
end


if isequal(sets.comparisonMethod,'mass_pcc')
    pvalfun = @(x,l1,l2) 1-beta_ev_cdf(x,alphaNu*l1,1,alphaN*2*l2,1);
    import Zeromodel.beta_ev_cdf; % pcc
    for i =1:length(comparisonStruct)
        if iscell(theoryStruct)
            comparisonStruct{i}.pval = pvalfun(comparisonStruct{i}.maxcoef(1),comparisonStruct{i}.lengthMatch,theoryStruct{comparisonStruct{i}.idx}.length);
        else
            comparisonStruct{i}.pval = pvalfun(comparisonStruct{i}.maxcoef(1),comparisonStruct{i}.lengthMatch,theoryStruct(comparisonStruct{i}.idx).length);
        end
    end
else % should run mp-pvalue or Stouffer. This has old p-value params / less accurate
    alphaNu = 0.001; %should be rather insensitive to specific value, check figure S5-S6;
   % p-value for local
    import Zeromodel.beta_ev_cdf; % correct form?
%     pvalfun = @(x,l1,l2,nuF,w) 1-beta_ev_cdf(x,nuF*w,1,nuF*2*(max(l1,l2)-w+1),0);
    pvalfun = @(x,l1,l2,nuF,w) 1-beta_ev_cdf(x,nuF*w,1,max(0.2*2*(max(l1,l2)-w),alphaNu*2*(l1-w+1)*(l2-w+1)),0);

    for i =1:length(comparisonStruct)
        if iscell(theoryStruct)
            comparisonStruct{i}.pval = pvalfun(comparisonStruct{i}.maxcoef(1),comparisonStruct{i}.bestLength,theoryStruct{comparisonStruct{i}.idx}.length,alphaNu, comparisonStruct{i}.lengthMatch);
        else
            comparisonStruct{i}.pval = pvalfun(comparisonStruct{i}.maxcoef(1),comparisonStruct{i}.bestLength,theoryStruct(comparisonStruct{i}.idx).length,alphaNu,comparisonStruct{i}.lengthMatch);
        end
    end
end

% % % 
% sets.timeFramesNr = nan;
% sets.displayResults=1
% sets.userDefinedSeqCushion = 0;
% sets.genConsensus = 0;
% import CBT.Hca.UI.get_display_results;
% [res] = get_display_results(barcodeGen,[], comparisonStruct, theoryStruct, sets);
% % % 
% % 
% % % % %% MP
% % numWorkers = 30;
% % MIN_OVERLAP_PIXELS_BL = 200;
% % minLen_BL = 200;
% % % minLen = 500; %150kb? or less? depends on application // if GUI, user selects thi
% % timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% % out=strcat('output',timestamp);
% % mkdir(out);
% % barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
% % [a,b] = find(barLens>minLen_BL);
% % % find(b==44)
% % import Core.compare_mp_all
% % ix = 1;
% % [mpI1individual,mp1individual,maxMP,stridxindividual,compStr] = ...
% %     compare_mp_all(theoryStruct,barcodeGen,minLen_BL,ix, timestamp,sF,MIN_OVERLAP_PIXELS_BL,numWorkers);
% % 
% % lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));
% % % lengthBorders = lengthBorders(ix);
% % % fig1 = figure;
% % bpPx = nmpx/nmbp;
% % import CBT.Hca.UI.Helper.plot_best_pos;
% % f= plot_best_pos([], compStr, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
% % % % % 
% % % % saveas(f,'figs/fig10.png')
% % % 
% % % % plot_best_pos([], comparisonStruct, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
% % 
% % 
% % %% for individual barcode
% % ix = 1;
% % compStr{ix}.allposB
% % 
% % strBestPos = cellfun(@(x) std(x.allposB(1:50)),compStr)
% % 
% % %where  strBestPos > threshold, possible false positive matches
% % figure,plot(strBestPos)


end
% % 
