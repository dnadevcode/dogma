function [compI] = generate_sf_struct(rezMaxMP,sets)

    allCoefs = cellfun(@(x) x{1}, rezMaxMP,'un',false);
    matAllCoefs =  cat(3, allCoefs{:});
  
    compI = struct('maxcoef',[],'pos',[],'or',[],'secondPos',[],'bestBarStretch',[],'bestlength',[]);

    %% Now find the best coefficient from matAllCoefs (using a cascading or whatever scheme)
    for barid =1:size(matAllCoefs,1)
        [singleCoef , singlePos ] =  max(matAllCoefs(barid,:,:),[],2);
        pos  = squeeze(singlePos)';
        compI.maxcoef(barid) =  squeeze(singleCoef);
        compI.or(barid) = arrayfun(@(x,y) rezMaxMP{x}{3}(barid,y), 1:length(rezMaxMP),pos);
        compI.pos(barid) = arrayfun(@(x,y) rezMaxMP{x}{2}(barid,y), 1:length(rezMaxMP),pos);
        compI.secondPos(barid) = arrayfun(@(x,y) rezMaxMP{x}{4}(barid,y), 1:length(rezMaxMP),pos);
        compI.bestlength(barid) = arrayfun(@(x,y) rezMaxMP{x}{5}(barid,y), 1:length(rezMaxMP),pos);
        compI.bestBarStretch(barid) =  sets.theory.stretchFactors(pos);
    end

end

