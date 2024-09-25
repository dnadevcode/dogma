function [] = quick_plot(idx,refBarcode,bars,origPos, origFlip,origStr)
% version of quick_visual_plot when theoryStruct needs to be loaded

% imresize(barcodeGen{idx}.rawBarcode,'Scale' ,[0 comparisonStruct{idx}.bestBarStretch])
curBar = imresize(bars{idx}.rawBarcode(logical(bars{idx}.rawBitmask)),'Scale' ,[1 1/origStr]) ;

if origFlip ~= 1
    curBar = fliplr(curBar);
end

thr = refBarcode;
thr = [thr nan(1,length(thr))];
% thrLambda = importdata(theoryStruct2{1}.filename);

% rawBg = barcodeGen{idx}.rawBg; %mode(cellfun(@(x) x.rawBg,barcodeGenLambda));
% meanlambda = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGenLambda));
% meanbar = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGen));
pos = find(bars{idx}.rawBitmask==1,1,'first');
figure,plot(pos+[origPos:origPos+length(curBar)-1],zscore(curBar))
hold on
plot(origPos:origPos+length(curBar)-1,zscore(thr(origPos:origPos+length(curBar)-1)))
legend({'experiment','thry'})

zscore(curBar,1)*zscore(thr(origPos:origPos+length(curBar)-1),1)'/length(curBar)
end

