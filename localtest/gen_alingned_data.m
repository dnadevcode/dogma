function [kymoStructs] = gen_alingned_data(kymoStructs)

sets.minOverlap = 200;
sets.maxShift = 20;
sets.skipPreAlign = 0;
sets.detPeaks = 1;

import OptMap.KymoAlignment.SPAlign.spalign;

% alignment. Look for features, so does not matter if we found the edges 
for idFold = 1:length(kymoStructs)    
%     idFold
    for idKym=1:length(kymoStructs{idFold})  % todo: should get directly from dbmstruct
        try
          [kymoStructs{idFold}{idKym}.alignedKymo,kymoStructs{idFold}{idKym}.alignedMask,cutKymo,cutMaskF,newMask] = ...
        spalign(double(kymoStructs{idFold}{idKym}.unalignedKymo),kymoStructs{idFold}{idKym}.unalignedBitmask,sets.minOverlap,sets.maxShift,sets.skipPreAlign, sets.detPeaks);
        catch
            kymoStructs{idFold}{idKym}.alignedKymo = [];
            kymoStructs{idFold}{idKym}.alignedMask  = [];
            newMask = [];

        end
          try
            kymoStructs{idFold}{idKym}.leftEdgeIdxs = arrayfun(@(frameNum) find(kymoStructs{idFold}{idKym}.alignedMask(frameNum, :), 1, 'first'), 1:size(kymoStructs{idFold}{idKym}.alignedMask,1));
            kymoStructs{idFold}{idKym}.rightEdgeIdxs = arrayfun(@(frameNum) find(kymoStructs{idFold}{idKym}.alignedMask(frameNum, :), 1, 'last'), 1:size(kymoStructs{idFold}{idKym}.alignedMask,1));
        catch
            kymoStructs{idFold}{idKym}.leftEdgeIdxs = [];
            kymoStructs{idFold}{idKym}.rightEdgeIdxs = [];
        end     
%         kymoStructs{idFold}{idKym}.name = 'somename';
        % this can be moved to be calculated separately
        kymoStructs{idFold}{idKym}.len = sum(newMask');
        if length( kymoStructs{idFold}{idKym}.len)>2
            kymoStructs{idFold}{idKym}.params =  polyfit(1:length( kymoStructs{idFold}{idKym}.len), kymoStructs{idFold}{idKym}.len,1);
        else
            kymoStructs{idFold}{idKym}.params = nan;
        end
    end

end

end

