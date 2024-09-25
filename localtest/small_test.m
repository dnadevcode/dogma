compStr{find(b==yId)}


thr = importdata(theoryStruct{1}.filename);

figure,plot(thr)

bar = barcodeGen{yId}.rawBarcode(barcodeGen{yId}.rawBitmask);


numWorkers = 30;
MIN_OVERLAP_PIXELS_BL = 200;
minLen_BL = 500;
% minLen = 500; %150kb? or less? depends on application // if GUI, user selects thi
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
[a,b] = find(barLens>minLen_BL);
find(b==44)
import Core.compare_mp_all
newB =barcodeGen(yId);
newB{1}.rawBarcode = fliplr(newB{1}.rawBarcode);
newB{1}.rawBitmask = fliplr(newB{1}.rawBitmask);

ix = 1;
[mpI1individual,mp1individual,maxMP,stridxindividual,compStr] = ...
    compare_mp_all(theoryStruct,newB,minLen_BL,ix, timestamp,1.01,MIN_OVERLAP_PIXELS_BL,numWorkers);
figure,plot(mpI1individual{1})
figure,plot(mp1individual{1})

%% PCC 
parpool('local')

parfor ix=1:12;
    numWorkers = 1;
    tic
    import Core.calc_overlap;
    [mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,{names2{ix}},{namesBar{ix}});
    toc
end
