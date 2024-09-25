meanNfrag = [];

fragSizes =  20:10:50;

for NUM_RAND_FRAGMENTS = fragSizes; % number of random fragments

    import Nullmodel.gen_rand;
    [barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr] = gen_rand(NUM_RAND_FRAGMENTS);

    lens = cellfun(@(x) x.lengthMatch,synthStr);

    probI = zeros(1,length(lens));

    minL = 150;

    for i=1:length(lens)
        probI(i) = (1-NUM_RAND_FRAGMENTS/theoryStruct{1}.length)^(lens(i)-minL-1); %minL-1 conseq
    end

    meanNfrag = [meanNfrag sum(probI)];
end

figure,plot(fragSizes,meanNfrag)