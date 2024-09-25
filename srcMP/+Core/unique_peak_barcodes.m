function [PKS, LOCS, pksUnique, pksUniquePos, bar1num, resfac2idx, or2sign] = unique_peak_barcodes(mp, mpI, baridx2, stridx, k)
% unique_peak_barcodes

%   Args:
%       mp - matrix profile
%       mpI - matrix profile index
%       baridx2 - bar1 nr on timeseries1 (non-stretched) 
%       stridx - re-scaling factor on bar2 (stretched
%       k - index of which comparison are we running
%

%   Returns:
%       PKS - mp values sorted
%       LOCS - locations of sorted PKS (mp values) along bar1
%       pksUnique - unique nr of bar1
%       pksUniquePos - indexes of the unique nr on bar1num
%       bar1num - specific barcode nr. for each location
%       resfac2idx - rescale factor of bar2    
%       or2sign - orientation sign of bar2



% [PKS,LOCS] = findpeaks(mp,'SortStr','descend','MinPeakHeight',thresh);
% [PKS,LOCS] = findpeaks(mp,'SortStr','descend','MinPeakHeight',thresh);

    % Sort MP values
    [PKS,LOCS] = sort(mp,'descend');
    
    % Remove any nan's (happens in between barcodes)
    LOCS(isnan(PKS)) = [];
    PKS(isnan(PKS)) = [];


    % which rescaling factor
    resfac2idx = stridx{k}(mpI(LOCS)+1);
    %
    or2sign = sign(resfac2idx);

    % which barcode do we match bar 1?
    bar1num = baridx2{k}(LOCS);
    [pksUnique,pksUniquePos] = unique(bar1num,'stable');


end

