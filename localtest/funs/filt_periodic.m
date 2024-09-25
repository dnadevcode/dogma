function [barcodeGenNonCirc, nonper] = filt_periodic(barcodeGen,timestamp,sF,ix)
    % filter out barcodes that match well to a periodic feature

barPeriodic  = barcodeGen{ix}.rawBarcode(barcodeGen{ix}.rawBitmask);

i=1;theoryStruct=[];
theoryStruct{i}.filename = fullfile(strcat('output',timestamp),'theory_barcode.txt');
fileID = fopen(theoryStruct{i}.filename,'w');
fprintf(fileID,strcat([' %2.' num2str(14) 'f ']), barPeriodic);
fclose(fileID);
theoryStruct{i}.meanBpExt_nm = 0.225;
theoryStruct{i}.psfSigmaWidth_nm = 300;
theoryStruct{i}.length = length(barPeriodic);
theoryStruct{i}.isLinearTF = 0;
theoryStruct{i}.name = 'Synthetic theory';


% figure,plot(barPeriodic)


% % %% MP
numWorkers = 30;
MIN_OVERLAP_PIXELS_BL = 100;
minLen_BL = 100;
% minLen = 500; %150kb? or less? depends on application // if GUI, user selects thi
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
[a,b] = find(barLens>minLen_BL);
% find(b==44)
import Core.compare_mp_all
ix = 1;
[mpI1individual,mp1individual,maxMP,stridxindividual,compStr] = ...
    compare_mp_all(theoryStruct,barcodeGen,minLen_BL,ix, timestamp,sF,MIN_OVERLAP_PIXELS_BL,numWorkers);

lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));

scores = cellfun(@(x) x.maxcoef, compStr);

scoreThresh = multithresh(scores);

barcodeGenNonCirc = barcodeGen(scores<scoreThresh);

nonper = scores<scoreThresh;
% 
end

