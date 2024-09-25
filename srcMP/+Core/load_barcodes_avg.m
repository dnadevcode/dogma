function [barcodeGen,passingKymos,kymoStructs,kS] = load_barcodes_avg(sets, timeFrameIdx,numBars)

if nargin < 2
    timeFrameIdx = sets.timeFrameIdx;
end

if nargin < 3
    numBars = 5;
end

passingKymos=[];

try mkdir(sets.output.matDirpath);catch; end
% add kymographs
import TempStuff.add_kymographs_fun;%CBT.Hca.Import.add_kymographs_fun;
[kymoStructs] = add_kymographs_fun(sets);

j=timeFrameIdx;
kS = cellfun(@(y) arrayfun(@(x) y.unalignedKymo(x:x+j,:),1:(timeFrameIdx+1):min(numBars,(size(y.unalignedKymo,1)-j)),'un',false),kymoStructs,'un',false);

lenBars = cellfun(@(x) length(x),kS);
idx = cell2mat(arrayfun(@(x) repmat(x,1,lenBars(x)),1:length(lenBars),'un',false));


kS =  [kS{:}];

for i=1:length(kS)
    kymoS{i}.unalignedKymo = kS{i};
    kymoS{i}.name = idx(i);
end

% align kymos - should not take any time if it's just single frame
% sets.alignMethod = 0;
import TempStuff.align_kymos; %CBT.Hca.Core.align_kymos;
[kymoStructsAligned] = align_kymos(sets,kymoS);
      
% ge

% 
% kymoStructs{1}.unalignedKymo(1:2,:)
% % % now seperate kymostructs based on timeFrameIdx (i.e. how many timeframes
% % % in a single one)
% % % 
% % % buffer is a matlab function that divides the data into 
% % % pieces with the length of k and with overlaps of m-1 
% % % samples and make a matrix s from the data vector. 
% % s=buffer(kymoStructs{1}.unalignedKymo,2,1); 
%             
%             
% % cut into INDEPENDENT slices. Using s
% kymoStructs{1}.unalignedKymo

% 
% %  put the kymographs into the structure
% import CBT.Hca.Core.edit_kymographs_fun;
% [kymoStructs, passingKymos] = edit_kymographs_fun(kymoStructs,sets.timeFramesNr,timeFrameIdx);
% 
% import generate.kymo_to_multi_bar;
% kymoStructs = kymo_to_multi_bar(kymoStructs);
% now convert this to many kymo of single frame

% align kymos - should not take any time if it's just single frame
% sets.alignMethod = 0;
% import CBT.Hca.Core.align_kymos;
% [kymoStructsAligned] = align_kymos(sets,kymoStructs);
%       
% generate barcodes / could generate multidimensional barcodes (i.e. each
% timeframe a different dimension
import TempStuff.gen_barcodes%CBT.Hca.Core.gen_barcodes;
barcodeGen =  TempStuff.gen_barcodes(kymoStructsAligned, sets);
 

for i=1:length(barcodeGen)
    barcodeGen{i}.kymo = kymoStructsAligned{i}.alignedKymo;
    barcodeGen{i}.bpsPerPx_original = -1;
end

end

