function [tests] = approx_main_kymo_molecule_edges_test()

tests = functiontests(localfunctions);

end


function test1Case(testCase)


       spaceList = {'Otsu','Double tanh','Error function','Zscore'}; 
            
%             [answer, tf] = listdlg('ListString', spaceList,...
%             'SelectionMode', 'Single', 'PromptString', 'Select item', 'Initialvalue', 1,'Name', 'Make choice');

    answer = 1;
    sets.skipDoubleTanhAdjustment = 1;

    if answer ~=1
        sets.edgeDetectionSettings.method = spaceList{answer}; 
    else
         sets.edgeDetectionSettings.method = 'Otsu';
        import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
        sets.edgeDetectionSettings = get_default_edge_detection_settings(sets.skipDoubleTanhAdjustment);
        % user canceled or closed dialog
    end

    kymoMatrix = zeros(10,200);
    kymoMatrix(:,20:180) = 1;

    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    [leftEdgeIdxs,rightEdgeIdxs,alignedMask] = approx_main_kymo_molecule_edges(kymoMatrix,  sets.edgeDetectionSettings);       
    
%     abs(mean(leftEdgeIdxs)-20)<=3
%     abs(mean(rightEdgeIdxs)-180)<=3
    verifyEqual(testCase,abs(mean(leftEdgeIdxs)-20)<=3,true);
    verifyEqual(testCase,abs(mean(rightEdgeIdxs)-180)<=3,true);

end



function test2Case(testCase)


       spaceList = {'Otsu','Double tanh','Error function','Zscore'}; 
            
%             [answer, tf] = listdlg('ListString', spaceList,...
%             'SelectionMode', 'Single', 'PromptString', 'Select item', 'Initialvalue', 1,'Name', 'Make choice');

    answer = 3;
    sets.skipDoubleTanhAdjustment = 1;

    if answer ~=1
        sets.edgeDetectionSettings.method = spaceList{answer}; 
    else
         sets.edgeDetectionSettings.method = 'Otsu';
        import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
        sets.edgeDetectionSettings = get_default_edge_detection_settings(sets.skipDoubleTanhAdjustment);
        % user canceled or closed dialog
    end

    kymoMatrix = zeros(10,200);
    kymoMatrix(:,20:180) = 1;

    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    [leftEdgeIdxs,rightEdgeIdxs,alignedMask] = approx_main_kymo_molecule_edges(kymoMatrix,  sets.edgeDetectionSettings);       
    
%     abs(mean(leftEdgeIdxs)-20)<=3
%     abs(mean(rightEdgeIdxs)-180)<=3
    verifyEqual(testCase,abs(mean(leftEdgeIdxs)-20)<=3,true);
    verifyEqual(testCase,abs(mean(rightEdgeIdxs)-180)<=3,true);

end



function test3Case(testCase)


       spaceList = {'Otsu','Double tanh','Error function','Zscore'}; 
            
%             [answer, tf] = listdlg('ListString', spaceList,...
%             'SelectionMode', 'Single', 'PromptString', 'Select item', 'Initialvalue', 1,'Name', 'Make choice');

    answer = 4;
    sets.skipDoubleTanhAdjustment = 1;

    if answer ~=1
        sets.edgeDetectionSettings.method = spaceList{answer}; 
    else
         sets.edgeDetectionSettings.method = 'Otsu';
        import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
        sets.edgeDetectionSettings = get_default_edge_detection_settings(sets.skipDoubleTanhAdjustment);
        % user canceled or closed dialog
    end

    kymoMatrix = zeros(10,200);
    kymoMatrix(:,20:180) = 1;

    import OptMap.MoleculeDetection.EdgeDetection.approx_main_kymo_molecule_edges;
    [leftEdgeIdxs,rightEdgeIdxs,alignedMask] = approx_main_kymo_molecule_edges(kymoMatrix,  sets.edgeDetectionSettings);       
    
%     abs(mean(leftEdgeIdxs)-20)<=3
%     abs(mean(rightEdgeIdxs)-180)<=3
    verifyEqual(testCase,abs(mean(leftEdgeIdxs)-20)<=3,true);
    verifyEqual(testCase,abs(mean(rightEdgeIdxs)-180)<=3,true);

end

