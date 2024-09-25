function [tests] = plot_best_bar_test()

tests = functiontests(localfunctions);

end


function test1Case(testCase)


    f = figure('Visible','off');
    ax = subplot(1,1,1);
    hold on
    import CBT.Hca.UI.Helper.plot_best_bar;

    barcodeGen{1}.rawBarcode = ones(1,300);
    barcodeGen{1}.rawBarcode(150:155) = 4;
    barcodeGen{1}.rawBitmask = ones(1,300);

    comparisonStruct{1}.idx =1;
    comparisonStruct{1}.pos =1;
    comparisonStruct{1}.or =1;
    comparisonStruct{1}.bestBarStretch = 1;
    
    theoryStruct(1).name = 'example';
    theoryStruct(1).rawBarcode =  1.1*ones(1,400);
        theoryStruct(1).rawBarcode(150:160) = 2;
     theoryStruct(1).rawBitmask = [];
      theoryStruct(1).length = 400;
     maxcoef = 0.9;

    userDefinedSeqCushion = 0;

    plot_best_bar(ax, barcodeGen, [], comparisonStruct,theoryStruct, maxcoef,userDefinedSeqCushion )
%     plot_best_bar
%     verifyEqual(testCase,abs(mean(rightEdgeIdxs)-180)<=3,true);




end
