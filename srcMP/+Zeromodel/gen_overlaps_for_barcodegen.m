function [] = gen_overlaps_for_barcodegen(sF,inputArg2)


sF = 0.8:0.01:1.2;

import Nullmodel.gen_random;
[~,barcodeGen] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,0);
%

minOverlap = 300;
tic
import Core.calc_overlap_pcc;
[overlapStruct] = calc_overlap_pcc(barcodeGen, sF,minOverlap);
toc

end

