function [overlapStruct2] = compare_mp_multi_theories_fast(barStruct, theoryStruct, sF, w, numWorkers)
% compare_mp_multi_theories_fast faster version where all theories are
% concatenated
%
%   Args:
%
%   Returns:
%       compStrOut, mp1individualOut
%

if nargin < 5
    numWorkers = 30;
end

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);

% lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);

if ~isstruct(barStruct)
    barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barStruct,'un',false);...
    cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
end

foldSynth = 'barcoli';


[names2, baridx2] = Core.save_long_bar_txt_skip(theoryStruct,inf,foldSynth);

baridx2{1} = baridx2{1}(1:end-w+2);

% save re-scaled exp. Single
[namesBar, stridx, ~] = Core.save_bars_rescaled_txt_stack(barStruct,sF,foldSynth);


stridx = stridx(1:end-w+1);
% baridx = baridx(1:end-w+1);

import Core.calc_overlap;

import Core.mp_to_struct_multi;

for ix=1:length(namesBar)

    [mpI1,mp1,maxMP] = calc_overlap([],timestamp, [], w, numWorkers,names2,namesBar);

    % now we need to convert from mp1 and mpI1 to usable table
    [overlapStruct2] = mp_to_struct_multi(length(theoryStruct), mp1,mpI1,baridx2,stridx,w,sF);
  

end

plotexample = 0;
if plotexample
% PLOT. First calculate full overlap
ix = 80;
theoryStruct(ix).rawBitmask = logical(ones(1,length(theoryStruct(ix).rawBarcode ) ));
import Core.get_full_overlap_score;
[overlapStruct2(1,ix).fullscore,overlapStruct2(1,ix).overlaplen, overlapStruct2(1,ix).lenB , overlapStruct2(1,ix).lenA,overlapStruct2(1,ix).partialScore,...
     overlapStruct2(1,ix).partialLength] = get_full_overlap_score(overlapStruct2(1,ix).pA,overlapStruct2(1,ix).pB,...
     overlapStruct2(1,ix).bestBarStretch, overlapStruct2(1,ix).or,[barStruct(1) theoryStruct(ix)],w);


%f= figure;
% tiledlayout(2,1)
import Plot.pair_evaluation_with_ground_truth_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_with_ground_truth_plot([barStruct(1) theoryStruct(ix)], [overlapStruct2(1,1) overlapStruct2(1,ix)],1,2,[],10000);
end



end

