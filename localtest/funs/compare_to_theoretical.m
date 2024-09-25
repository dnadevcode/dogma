function [comparisonStruct] = compare_to_theoretical(barcodeGen,theoryStruct,sF,sets)

%% PCC
% barcodeGen = barcodeGen1;
if isempty(sF)
    sF = 0.8:0.01:1.2;
end

sets.filterSettings.filter = 0;

rezMax=[];bestBarStretch=[];bestLength=[];
for i=1:length(theoryStruct)
    tic
    import CBT.Hca.Core.Comparison.on_compare;
    [rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGen,...
        theoryStruct{i},sets.comparisonMethod, sF,[],50,[],sets.filterSettings);
    toc
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
