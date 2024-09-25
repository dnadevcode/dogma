function [tests] = run_pipeline_scores_test()
tests = functiontests(localfunctions);
end


function test_run_pipeline_scores_1Case(testCase)

    tempFold1 = tempname;
    fullFold = fullfile(tempFold1,'test1','211213_Sample358-3-st2_647.81bpPERpx_0.169nmPERbp','kymos');
    [~,~] = mkdir(fullFold);
    tempMat = zeros(1,400);
    tempMat([100:299]) = 10+normrnd(0,1,1,200);
    tempMat = repmat(tempMat,20,1);
    tempMat = imgaussfilt(tempMat,[1 3]);
    imwrite(uint16((2^16-1)*(tempMat-min(tempMat(:)))/(max(tempMat(:))-min(tempMat(:)))),fullfile(fullFold,'kymo.tif'))

%     import Helper.load_kymo_data;
%     [kymoStructs,barGen] = load_kymo_data(tempFold1,1,1,1,1);


    import Helper.get_all_folders;
    [barN, twoList] = get_all_folders(tempFold1);

  
    w = [0 250];
    sF = 1;
    timeFramesNr = 0;

    %%theory       
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
    ligandLength = 4;

    import Helper.get_theory_twostate;
    [theorySeq,theoryStruct] = get_theory_twostate(fastaName,gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY, cN,kY,ligandLength);

    theoryGen = [];
    theoryGen.theoryBarcodes = {theoryStruct{1}.rawBarcode};
    theoryGen.theoryBitmasks = {theoryStruct{1}.rawBitmask};
    theoryGen.theoryNames = {theoryStruct{1}.name};
    theoryGen.theoryIdx =  {1};
%     theoryGen.sets = hcatheorySets;
    theoryGen.sets.meanBpExt_nm = 0.3;
    theoryGen.sets.pixelWidth_nm = nmpx;
    theoryGen.sets.psfSigmaWidth_nm = psf;
     theoryGen.sets.isLinearTF = ~isC;

    matFilename = strcat(['theoryGen_', num2str(0.3) '_' num2str(nmpx) '_' num2str(psf) '_' num2str(~isC) '_session.mat']);
    matFilepath = fullfile(tempFold1, matFilename);
    
    save(matFilepath,'theoryGen')
%     thryFiles 
    thryFiles = dir(fullfile(tempFold1,'*.mat'));
    
    [t] = run_pipeline_scores(tempFold1, twoList(1,:), 1, [0], sF, timeFramesNr, thryFiles);


    [~,~] =  rmdir(tempFold1,'s'); % remove temporary folder


%%
% 
%     thryFiles = dir(fullfile(tempFold1,'*.mat'));
% 
%     c = parcluster;
%     c.AdditionalProperties.AccountName = 'naiss2024-22-957';
%     c.AdditionalProperties.WallTime = '6:00:00';
% 
% 
%     sF = 0.8:0.025:1.2;
%     w= 0 ; 
%     numW = 60;
% 
%     % w = [ 500:50:1000];
%     idd = 7;
%     batch(c,@run_pipeline_scores,1,{dirName,[twoList(idd,:)], 1, w, sF, thryFiles},'Pool',numW);
% 
% 





end

