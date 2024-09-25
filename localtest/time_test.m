
    
vals = 5:10:50;
time = zeros(1,length(vals));
for numB = 1:length(vals)
    bars = bgAll(cellfun(@(x) sum(x.rawBitmask),bgAll)>sets.minOverlap); % only for those larger than min overlap
    bars = bars(1:vals(numB));


    tic

   
    [oS] = calc_overlap_mp(bars,sF, sets.minOverlap,timestamp);
time(numB) = toc

end
