function [barcodeGen,passingKymos,kymoStructs] = load_barcodes(sets, timeFrameIdx)

if nargin < 2
    timeFrameIdx = sets.timeFrameIdx;
end


try mkdir(sets.output.matDirpath);catch; end
% add kymographs
import CBT.Hca.Import.add_kymographs_fun;
[kymoStructs] = add_kymographs_fun(sets);

%  put the kymographs into the structure
import CBT.Hca.Core.edit_kymographs_fun;
[kymoStructs, passingKymos] = edit_kymographs_fun(kymoStructs,sets.timeFramesNr,timeFrameIdx);
% 
% import generate.kymo_to_multi_bar;
% kymoStructs = kymo_to_multi_bar(kymoStructs);
% now convert this to many kymo of single frame

% align kymos - should not take any time if it's just single frame
% sets.alignMethod = 0;
import CBT.Hca.Core.align_kymos;
[kymoStructsAligned] = align_kymos(sets,kymoStructs);
      
% generate barcodes / could generate multidimensional barcodes (i.e. each
% timeframe a different dimension
import CBT.Hca.Core.gen_barcodes;
barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructsAligned, sets);
 

for i=1:length(barcodeGen)
    barcodeGen{i}.kymo = kymoStructsAligned{i}.alignedKymo;
end

end

