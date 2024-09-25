nT = 150;
tpList =[];
for ii=1:nT
    idxPair = ii;
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    px1 = comparisonStruct{xId}.pos(1):comparisonStruct{xId}.pos(1)+comparisonStruct{xId}.lengthMatch-1;
    px2 = comparisonStruct{yId}.pos(1):comparisonStruct{yId}.pos(1)+comparisonStruct{yId}.lengthMatch-1;
    if length(overLap)>=100
%         tp(ii)=1;
                tpList = [tpList xId yId];

    end
end

% barsAlignment = unique(tpList);

  import CBT.Hca.UI.Helper.plot_best_pos;
%     f = plot_best_pos([], compStr([find(b==xId) find(b==yId)]), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%  saveas(f,'figs/fig12.png')
% 
    plot_best_pos([], comparisonStruct(barsAlignment), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%  