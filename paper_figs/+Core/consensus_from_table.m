function [outConsensus] = consensus_from_table(tableS,barcodeGen)


   %% plot multi theories with some nan's in between

   N = length(barcodeGen);

   lengthsbarCur = max(tableS(:,2))-min(tableS(:,1))+1;

   tableS(:,1:2) = tableS(:,1:2)-min(tableS(:,1))+1;

        outConsensus = nan(size(tableS,1),lengthsbarCur);
%     outConsensus(1,:) = zscore(theorBar);
%     lenT = length(theorBar);
    
        for k=1:N
            expBar = barcodeGen{k}.rawBarcode;
            expBit = barcodeGen{k}.rawBitmask;

            expBar = expBar(expBit);

            % interpolate to the length which gave best CC value
            expBar = imresize(expBar,'Scale' ,[1 tableS(k,4)]); % ?
            expBit = imresize(expBit,'Scale' ,[1 tableS(k,4)]);
            expBar(~expBit)= nan;

            if tableS(k,3)==2
                expBar = fliplr(expBar);
                expBit = fliplr(expBit);
            end
            expBar = expBar(expBit);

            outConsensus(k,tableS(k,1):tableS(k,1)+length(expBar)-1) = (expBar-nanmean(expBar))/nanstd(expBar,1)+5;
        end
% outConsensu


end

