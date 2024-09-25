function [] = quick_visual_plot(idx,thrIdx,bars,rezMax,bestBarStretch,theoryStruct)

% imresize(barcodeGen{idx}.rawBarcode,'Scale' ,[0 comparisonStruct{idx}.bestBarStretch])
curBar = imresize(bars{idx}.rawBarcode(bars{idx}.rawBitmask),'Scale' ,[1 bestBarStretch{thrIdx}(idx)]) ;

if rezMax{thrIdx}{idx}.or(1)==2
    curBar = fliplr(curBar);
end

thr = theoryStruct(thrIdx).rawBarcode;
thr = [thr nan(1,length(thr))];
% thrLambda = importdata(theoryStruct2{1}.filename);

% rawBg = barcodeGen{idx}.rawBg; %mode(cellfun(@(x) x.rawBg,barcodeGenLambda));
% meanlambda = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGenLambda));
% meanbar = mean(cellfun(@(x) mean(x.rawBarcode(x.rawBitmask)),barcodeGen));
pos = find(bars{idx}.rawBitmask==1,1,'first');
figure,plot(pos+[rezMax{thrIdx}{idx}.pos:rezMax{thrIdx}{idx}.pos+length(curBar)-1],zscore(curBar))
hold on
plot(rezMax{thrIdx}{idx}.pos:rezMax{thrIdx}{idx}.pos+length(curBar)-1,zscore(thr(rezMax{thrIdx}{idx}.pos:rezMax{thrIdx}{idx}.pos+length(curBar)-1)))
legend({'experiment','thry'})


end

