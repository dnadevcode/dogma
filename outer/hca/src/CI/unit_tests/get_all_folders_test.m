function tests = get_all_folders_test
    tests = functiontests(localfunctions);
end

% results = runtests('get_all_folders_test')


function testTestDatacase(testCase)

    % Imitates the structure of kymo folder
    tempFold1 = tempname;
    fullFold = fullfile(tempFold1,'test1','211213_Sample358-3-st2_647.81bpPERpx_0.169nmPERbp','kymos');
    [~,~] = mkdir(fullFold);


    import Helper.get_all_folders;
    [barN, twoList] = get_all_folders(tempFold1);

    [~,~] =  rmdir(tempFold1,'s'); % remove temporary folder

    verifyEqual(testCase,barN{1},0)
    verifyEqual(testCase,twoList, [1 1])


end
