function [tests] = hca_compare_distance_test()
    tests = functiontests(localfunctions);
end


function test1Case(testCase)



    import CBT.Hca.Core.Comparison.hca_compare_distance;

    sets.w = 0;
    sets.comparisonMethod = 'mpnan';

    N = 100;
    barGen{1}.rawBarcode = normrnd(0,1,1,N);
    barGen{1}.rawBitmask = ones(1,N);

    M = 200;

    theoryStruct(1).rawBarcode = normrnd(0,1,1,M);
    theoryStruct(1).rawBitmask = [];
    theoryStruct(1).name = 'test';
    theoryStruct(1).isLinearTF = 1;

  
    if sets.w == 0
        sets.comparisonMethod = 'mass_pcc';
    else
        sets.comparisonMethod = 'mpnan';
    end

    sets.theory.stretchFactors = [0.9 1 1.1];
    import Core.rescale_barcode_data;
    [barGen] = rescale_barcode_data(barGen,sets.theory.stretchFactors);


    [rezMax] = hca_compare_distance(barGen, theoryStruct, sets );



    verifyEqual(testCase,sum(isnan(rezMax{1}{1})),0);



end
