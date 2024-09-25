function [tests] = load_kymo_data_test()

    tests = functiontests(localfunctions);
    % results = runtests('load_kymo_data_test')

end

function testKymo1Case(testCase)

    tempFold1 = tempname;
    fullFold = fullfile(tempFold1,'test1','211213_Sample358-3-st2_647.81bpPERpx_0.169nmPERbp','kymos');
    [~,~] = mkdir(fullFold);
    tempMat = zeros(1,400);
    tempMat([100:299]) = 10+normrnd(0,1,1,200);
    tempMat = repmat(tempMat,20,1);
    tempMat = imgaussfilt(tempMat,[1 3]);
    imwrite(uint16((2^16-1)*(tempMat-min(tempMat(:)))/(max(tempMat(:))-min(tempMat(:)))),fullfile(fullFold,'kymo.tif'))

    import Helper.load_kymo_data;
    [kymoStructs,barGen] = load_kymo_data(tempFold1,1,1,1,1);


  
    [~,~] =  rmdir(tempFold1,'s'); % remove temporary folder

    verifyEqual(testCase,length(barGen),1);


end