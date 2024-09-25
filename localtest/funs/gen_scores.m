function [val] = gen_scores(fastaFile,nmbp,nmpx,psffac,nmbmfac,barcodeGen,sets,sF)


%     nmbp = 0.2;
% nmpx = 110; % 208?
nmbp = nmbp*psffac*nmbmfac;
nmpx = nmpx*psffac;

% idx = 1;
import Thry.gen_theoretical;
[theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmpx);

% sF = 1;
barcodeGenSub = barcodeGen;

rezMax=[];bestBarStretch=[];bestLength=[]; % todo: write this as overlapstructure?
for i=1:length(theoryStruct)
    tic
    import CBT.Hca.Core.Comparison.on_compare;
    [rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGenSub,...
        theoryStruct{i},sets.comparisonMethod, sF,[],50,[],sets.filterSettings);
    toc
end

comparisonStructAll = rezMax;
for i=1:length(comparisonStructAll)
    for j=1:length(bestBarStretch{i})
        comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
        comparisonStructAll{i}{j}.length = bestLength{i}(j);
    end
end
import CBT.Hca.Core.Comparison.combine_theory_results;
[comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);

val = mean(cellfun(@(x) x.maxcoef(1),comparisonStruct));


end

