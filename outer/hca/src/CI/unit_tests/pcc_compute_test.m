function tests = pcc_compute_test
    tests = functiontests(localfunctions);
end


% results = runtests('pcc_compute_test.m')

function testpcccase(testCase)
    % random test - run multiple times to catch edge cases
    rng('default')
    barcodeGen =[];
    barcodeGen{1}.rawBarcode = imgaussfilt(normrnd(0,1,1,1000),3);
    barcodeGen{1}.rawBitmask = ones(1,1000);
    barcodeGen{2}.rawBarcode = imgaussfilt(normrnd(0,1,1,1000),3);
    barcodeGen{2}.rawBitmask = ones(1,1000);
    theorBar = [barcodeGen{1}.rawBarcode imgaussfilt(normrnd(0,1,1,1000),3) ]; %imgaussfilt(normrnd(0,1,1,1000),3);
    theorBit = ones(1,1000);
    isLinearTF = 1;
    theoryStruct = bar_to_theory({theorBar theorBar}, {theorBit theorBit},isLinearTF);

    
    stretchFactors = 1;
    comparisonMethod = 'mass_pcc';
    w = 100;
    numPixelsAroundBestTheoryMask = 50;
    %      sets.theory.stretchFactors = 0.9:0.01:1.1;
    % import CBT.Hca.Core.Comparison.compute_distance;

    % [comparisonStruct,comparisonStructAll,rezMax,bestBarStretch,bestLength] = ...
    % compute_distance(barcodeGen,theoryStruct, setsThry, [],'mpnan' );
    import CBT.Hca.Core.Comparison.on_compare_barcode;
    [ rezMaxM,bestBarStretch,bestLength ] = on_compare_barcode(barcodeGen,theorBar, theorBit, isLinearTF, comparisonMethod,stretchFactors,w,numPixelsAroundBestTheoryMask);
    %       import CBT.Hca.Core.Comparison.combine_theory_results;
    %         [comparisonStruct] = combine_theory_results(theoryStruct, rezMaxM,bestBarStretch,bestLength);

    % for this should give rezMaxM as input
    rezMaxM{1}.idx = 1;
    rezMaxM{1}.bestBarStretch = bestBarStretch(1);
    rezMaxM{1}.bestLength = bestLength(1);
    rezMaxM{2}.idx = 1;
    rezMaxM{2}.bestBarStretch = bestBarStretch(2);
    rezMaxM{2}.bestLength = bestLength(2);
    setsThry = [];
    setsThry.displayResults = 1;
    setsThry.genConsensus = 0;
    setsThry.comparisonMethod = comparisonMethod;
    setsThry.theory.isLinearTF = isLinearTF;
    setsThry.timeFramesNr = 0;
    setsThry.plotallmatches = 0;
    setsThry.output.matDirpath = 'output';
    setsThry.timestamp =  datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    setsThry.w = 100;
    setsThry.userDefinedSeqCushion = 50;
    import CBT.Hca.UI.get_display_results;
    resStruct = get_display_results(barcodeGen,[], rezMaxM, theoryStruct, setsThry);

    verifyEqual(testCase,1,1)


% actSolution = quadraticSolver(1,-3,2);
% expSolution = [2 1];
end


function testpcccase2(testCase)
    % random test - run multiple times to catch edge cases
    rng('default')
    barcodeGen =[];
    barcodeGen{1}.rawBarcode = imgaussfilt(normrnd(0,1,1,1000),3);
    barcodeGen{1}.rawBitmask = ones(1,1000);
    barcodeGen{2}.rawBarcode = imgaussfilt(normrnd(0,1,1,1000),3);
    barcodeGen{2}.rawBitmask = ones(1,1000);
    theorBar = [barcodeGen{1}.rawBarcode imgaussfilt(normrnd(0,1,1,1000),3) ]; %imgaussfilt(normrnd(0,1,1,1000),3);
    theorBit = ones(1,1000);
    isLinearTF = 0;
    theoryStruct = bar_to_theory({theorBar theorBar}, {theorBit theorBit},isLinearTF);

    
    stretchFactors = 1;
    comparisonMethod = 'mass_pcc';
    w = 100;
    numPixelsAroundBestTheoryMask = 50;
    %      sets.theory.stretchFactors = 0.9:0.01:1.1;
    % import CBT.Hca.Core.Comparison.compute_distance;

    % [comparisonStruct,comparisonStructAll,rezMax,bestBarStretch,bestLength] = ...
    % compute_distance(barcodeGen,theoryStruct, setsThry, [],'mpnan' );
    import CBT.Hca.Core.Comparison.on_compare_barcode;
    [ rezMaxM,bestBarStretch,bestLength ] = on_compare_barcode(barcodeGen,theorBar, theorBit, isLinearTF, comparisonMethod,stretchFactors,w,numPixelsAroundBestTheoryMask);
    %       import CBT.Hca.Core.Comparison.combine_theory_results;
    %         [comparisonStruct] = combine_theory_results(theoryStruct, rezMaxM,bestBarStretch,bestLength);

    % for this should give rezMaxM as input
    rezMaxM{1}.idx = 1;
    rezMaxM{1}.bestBarStretch = bestBarStretch(1);
    rezMaxM{1}.bestLength = bestLength(1);
    rezMaxM{2}.idx = 1;
    rezMaxM{2}.bestBarStretch = bestBarStretch(2);
    rezMaxM{2}.bestLength = bestLength(2);
    setsThry = [];
    setsThry.displayResults = 1;
    setsThry.genConsensus = 0;
    setsThry.comparisonMethod = comparisonMethod;
    setsThry.theory.isLinearTF = isLinearTF;
    setsThry.timeFramesNr = 0;
    setsThry.plotallmatches = 0;
    setsThry.output.matDirpath = 'output';
    setsThry.timestamp =  datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    setsThry.w = 100;
    setsThry.userDefinedSeqCushion = 50;
    import CBT.Hca.UI.get_display_results;
    resStruct = get_display_results(barcodeGen,[], rezMaxM, theoryStruct, setsThry);

    verifyEqual(testCase,1,1)


% actSolution = quadraticSolver(1,-3,2);
% expSolution = [2 1];
end
