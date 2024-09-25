% triangle_test

% first we find the first element that is not significant
accepted = cellfun(@(x) x.accepted, barsetGen.currentData); % all the
firstNonAcc = find(~accepted,1,'first');

% based on this we can get the elements that are not accepted


f=figure,g=tiledlayout(5,2,'TileSpacing','loose')

% idxPair = 11;
[curSink,curSource] = ind2sub(size(oS),sortedIds(firstNonAcc));
import Plot.pair_evaluation_with_ground_truth_simple_plot; % todo: convert comparisonStructC to overlap for this plot (using something like synth_to_struct)
pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,curSink,curSource,[],[],g);


%% now triangle
f=figure,g=tiledlayout(5,2,'TileSpacing','loose')

pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,538,659,[],[],g);

pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,571,659,[],[],g);

pair_evaluation_with_ground_truth_simple_plot(barcodeGen, oS,538,571,[],[],g);


oS(538,571)

sub2ind(size(oS),571,659)
sortedIds(10)


[a b] = ind2sub(size(oS),sortedIds(10))

edge1 = [538 659]

edge2 = [571 659];

edge3 = [538 571];
nextElt = min(find(sortedIds == sub2ind(size(oS),edge3(1),edge3(2))),find(sortedIds == sub2ind(size(oS),edge3(2),edge3(1))))


barsetGen.currentData{nextElt}
% find(sortedIds == sub2ind(size(oS),edge3(2),edge3(1)))



find(barIslands{2}==538)
find(barIslands{2}==659)