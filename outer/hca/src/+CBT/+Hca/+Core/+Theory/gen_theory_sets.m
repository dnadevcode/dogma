function [sets] = gen_theory_sets(meanBpExt_nm,pixelWidth_nm,psfSigmaWidth_nm,isLinearTF,resultsDir)
    % import CBT.Hca.Core.Theory.gen_theory_sets

    %%
    sets.promptfortheory = 0;
    sets.promptforsavetheory = 0;
    sets.theoryGen.precision = 8;
    sets.skipBarcodeGenSettings = 1;
    sets.savetxts = 1;
    sets.theoryGen.meanBpExt_nm = meanBpExt_nm; % create bigger
    sets.theoryGen.concN = 60;
    sets.theoryGen.concY = 0.2;
    sets.theoryGen.concDNA = 0.2;
    sets.theoryGen.psfSigmaWidth_nm = psfSigmaWidth_nm;
    sets.theoryGen.pixelWidth_nm = pixelWidth_nm;
    sets.theoryGen.deltaCut = 3;
    sets.theoryGen.isLinearTF = isLinearTF;
    sets.theoryGen.widthSigmasFromMean= 4;
    sets.theoryGen.computeFreeConcentrations = 1;
    sets.theoryGen.k = 2^16;
    sets.theoryGen.m = 2^15;
    sets.theoryGen.method = 'literature';
    sets.model.pattern = '';
    sets.resultsDir = resultsDir;
    sets.lambda.fold = '';
    sets.lambda.name = 'sequence.fasta';
    sets.addgc = 0;
%     sets.output.matDirpath = 'output';
    import CBT.Hca.Core.Theory.choose_cb_model;
    [sets.model ] = choose_cb_model(sets.theoryGen.method,sets.model.pattern);

end

