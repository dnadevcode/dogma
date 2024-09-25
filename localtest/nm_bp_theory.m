function [outputArg1,outputArg2] = nm_bp_theory()

% can change: nmbp,nmpx, nmpsf( by changing ratio of other two


sets.timeFramesNr = nan;
sets.displayResults=1
sets.userDefinedSeqCushion = 0;
sets.genConsensus = 0;

sets.comparisonMethod = 'mass_pcc';
sets.filter = 0;
sets.filterSettings.filter = 0;

% psf = 300; %300 nm
val = [];

sF = 0.9:0.01:1.1;

fac = 0.8:0.01:1.2;
nmbmfac = 1;
psffac = 1;

nmbp = 0.22;
nmpx = 110; % 208?
val =[];
for idx = 1:length(fac);


 val(idx) = gen_scores(fastaFile, nmbp, nmpx, psffac, fac(idx), barcodeGen, sets, sF);

end

figure,plot(fac,val)
[a,b] = max(val)

nmFac = fac(b);


fac2 =  0.5:0.01:0.7;

valPSF =[];
for idx = 1:length(fac2);
    valPSF(idx) = gen_scores(fastaFile, nmbp, nmpx, fac2(idx), nmFac, barcodeGen, sets, sF);
end
figure,plot(fac2,valPSF)

[c,d]  = max(valPSF)
psfFac = fac2(d);

nmbpNew = nmbp*nmFac;
psfNew = 300/psfFac;



fac3 =  0.8:0.01:1.2;

valPX =[];
for idx = 1:length(fac2);
    valPX(idx) = gen_scores(fastaFile, nmbp, nmpx*fac3(idx), psfFac, nmFac, barcodeGen, sets, sF);
end

figure,plot(fac3(1:length(valPX)),valPX)
sets.timeFramesNr = nan;
sets.displayResults=1
sets.userDefinedSeqCushion = 0;
sets.genConsensus = 0;
import CBT.Hca.UI.get_display_results;
[res] = get_display_results(barcodeGenSub,[], comparisonStruct, theoryStruct, sets);

end

