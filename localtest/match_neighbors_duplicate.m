% sF = 0.8:.01:1.2;
sF = 1;
minOverlap = 100;
tic
oSneighbor = [];
for i=1:length(barcodeGen)-1
    import Core.calc_overlap_pcc_sort_m;
    [oSneighbor{i}] = calc_overlap_pcc_sort_m([barcodeGen(i:i+1)], sF,minOverlap);
end

potentialDuplicate = zeros(1,length(barcodeGen)-1);
scScore = nan(1,length(barcodeGen)-1);
for i=1:length(barcodeGen)-1
    potentialDuplicate(i) = oSneighbor{i}(1,2).sc > 0.7;
    if  potentialDuplicate(i)==1
        scScore(i) = oSneighbor{i}(1,2).sc;
    end
end

[a,b] = min(scScore);

potDubs = find(potentialDuplicate);
%%
% idx = b
% i = potDubs(idx);

i = b ;

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen(i:i+1),'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen(i:i+1),'un',false)]',{'rawBarcode','rawBitmask'},2);


import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, oSneighbor{i},1, 2,barStruct,nan);
%%
% average the two barcodes

pB = oSneighbor{i}(1,2).pA;
pA = oSneighbor{i}(1,2).pB;
lenA = length(barcodeGen{i}.rawBarcode);
lenB = length(barcodeGen{i+1}.rawBarcode);

stIdx = min(pA,pB);
pA  = pA-stIdx+1;
pB = pB-stIdx+1;
stopIdx =max(pA+lenA-1,pB+lenB-1);
twoBars = nan(2,stopIdx-stIdx+1);

barcodeGen{i}.rawBarcode(~barcodeGen{i}.rawBitmask) = nan;
barcodeGen{i+1}.rawBarcode(~barcodeGen{i+1}.rawBitmask) = nan;

twoBars(1,pA:pA+lenA-1)=barcodeGen{i}.rawBarcode;
twoBars(2,pB:pB+lenB-1)=barcodeGen{i+1}.rawBarcode;

% figure 
% plot(pA:pA+lenA-1,barcodeGen{i}.rawBarcode)
% hold on
% plot(pB:pB+lenB-1,barcodeGen{i+1}.rawBarcode)

figure,plot(twoBars')
newBar = nanmean(twoBars);
newBit = isnan(newBar);