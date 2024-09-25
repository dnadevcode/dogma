function [maxCoef] = calc_coefs(theoryStruct,barGen,windowWidths,sF,sets)


% windowWidths = 400:100:600;
sets.comparisonMethod = 'mpnan';

thryNames = {theoryStruct.name};
barcodeNames = cellfun(@(x) x.name,barGen,'un',false);

import CBT.Hca.Core.Comparison.hca_compare_distance;
import Core.export_coefs_local;

for wIdx = 1:length(windowWidths)
    sets.w = windowWidths(wIdx);
    display(['Running w = ',num2str(sets.w)]);
    % only local lengths for which all length re-scaled versions passes the
    % threshold
    passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),barGen)*sF(1) >= sets.w);

    if sets.w == 0
        sets.comparisonMethod = 'mass_pcc';
    else
        sets.comparisonMethod = 'mpnan';
    end

    [rezMaxMP] = hca_compare_distance(barGen(passingThreshBars), theoryStruct, sets );
    
    % Get info for all the coefficients. Save this as a separate matrix
    allCoefs = cellfun(@(x) x{1}, rezMaxMP,'un',false);
    matAllCoefs =  cat(3, allCoefs{:});
    save([sets.dirName, num2str(sets.w),'_', num2str(min(sF)),'sf_allcoefs.mat'],'matAllCoefs','passingThreshBars','sets','-v7.3');

    tic
    maxCoef = cell(1,size(matAllCoefs,1));
    maxOr = cell(1,size(matAllCoefs,1));
    maxPos = cell(1,size(matAllCoefs,1));
    maxSecondPos = cell(1,size(matAllCoefs,1));
    maxlen = cell(1,size(matAllCoefs,1));
    bestSF = cell(1,size(matAllCoefs,1));

%% Now find the best coefficient from matAllCoefs (using a cascading or whatever scheme)
    for barid =1:size(matAllCoefs,1)
        [singleCoef , singlePos ] =  max(matAllCoefs(barid,:,:),[],2);
        pos  = squeeze(singlePos)';
        maxCoef{barid} =  squeeze(singleCoef);
        maxOr{barid} = arrayfun(@(x,y) rezMaxMP{x}{3}(barid,y), 1:length(rezMaxMP),pos);
        maxPos{barid} = arrayfun(@(x,y) rezMaxMP{x}{2}(barid,y), 1:length(rezMaxMP),pos);
        maxSecondPos{barid} = arrayfun(@(x,y) rezMaxMP{x}{4}(barid,y), 1:length(rezMaxMP),pos);
        maxlen{barid} = arrayfun(@(x,y) rezMaxMP{x}{5}(barid,y), 1:length(rezMaxMP),pos);
        bestSF{barid} =  sets.theory.stretchFactors(pos);
    end
    toc
  
    export_coefs_local(thryNames,maxCoef,maxOr,maxPos,maxlen, bestSF, barcodeNames(passingThreshBars),[sets.dirName, '_MP_w=',num2str(sets.w),'_']);




end

end

