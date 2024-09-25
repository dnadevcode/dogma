function [PKS1,LOCS1,pksUnique1,pksUniquePos1,barcodePair1,rescaleFactorPair1,orPair1] =...
    calc_peak_locs(mp1,mpI1,baridx2,stridx,minPCC)

    %   Args:
    %       mp1,mpI1,baridx2,stridx,minPCC
    %
    %   Returns:
    %       PKS1,LOCS1,pksUnique1,pksUniquePos1,barcodePair1,rescaleFactorPair1,orPair1
    %

    NN=length(mpI1);
    PKS1=cell(1,NN);
    LOCS1=cell(1,NN);
    pksUnique1 = cell(1,NN);
    pksUniquePos1=cell(1,NN);
    barcodePair1=cell(1,NN);
    rescaleFactorPair1=cell(1,NN);
    orPair1=cell(1,NN);
    pvals =zeros(1,NN);
    for k=1:NN
%         k

        [PKS1{k},LOCS1{k},pksUnique1{k},pksUniquePos1{k},barcodePair1{k},rescaleFactorPair1{k},orPair1{k}] = ...
        Core.unique_peak_barcodes(mp1{k},mpI1{k}, minPCC,baridx2,stridx,k);
    end

end

