function [tests] = gen_simple_theory_px_test()

    tests = functiontests(localfunctions);
    % results = runtests('gen_simple_theory_px_test')

end

function test1Case(testCase)

    tempFold1 = tempname;
    [~,~] = mkdir(tempFold1);
        
    fastaName = fullfile(tempFold1,'tmp.fasta');
    fastawrite(fastaName,randseq(20000));


    sigma = 0.68; % scaling of AT-GCs
    gcSF  = 1;
    kY = 10; % binding constants, hard-coded for now
    kN = 30;
    psf = 370; % seemed to give best results
    cN = 6;
    cY = 0.04;
    nmpx = 110;
    isC = 1;
    pxSize = nmpx/0.3;



    import Helper.get_theory_twostate;
    [theorySeq,theoryStruct] = get_theory_twostate(fastaName,gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY, cN,kY);


    [~,~] =  rmdir(tempFold1,'s'); % remove temporary folder

    verifyEqual(testCase,length(theorySeq),54);


end
