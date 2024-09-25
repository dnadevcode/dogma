function [sets] = non_default_settings(sets)
    %
%     import CBT.Hca.Settings.non_default_settings

        
        % add filter settings if filtering is enabled
        if sets.filterSettings.filter
            % here we add settings
            prompt = {'Filter method','Filter size (px)'};
            title = 'Filtering settings';
            dims = [1 35];
            definput = {'0','2.3'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.filterSettings.filterMethod = str2double(answer{1}); % filter before or after stretching
            sets.filterSettings.filterSize = str2double(answer{2});      
        end

        if sets.genConsensus
           % here we add settings for this 
            prompt = {'Averaging method','Prompt for threshold', 'Consensus threshold'};
            title = 'Consensus settings';
            dims = [1 35];
            definput = {'bgmean','1','0.75'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.consensus.barcodeNormalization = answer{1};
            sets.consensus.promptForBarcodeClusterLimit =  str2double(answer{2}); 

            sets.consensus.threshold =  str2double(answer{3}); 
        end

        if sets.bitmaskSettings
            % here we add settings
            prompt = {'Camera (nm/px)','Point spread function (nm)', 'Delta cut'};
            title = 'Bitmasking settings';
            dims = [1 35];
            definput = {'130','300','3'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.bitmasking.prestretchPixelWidth_nm = str2double(answer{1}); % camera nm/px
            sets.bitmasking.psfSigmaWidth_nm = str2double(answer{2}); % psf in nanometers
            sets.bitmasking.deltaCut = str2double(answer{3}); % how many delta's should we take

            % number of untrusted pixels
            sets.bitmasking.untrustedPx = sets.bitmasking.deltaCut*sets.bitmasking.psfSigmaWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;
        end
    
        % edge settings. Add a new possibility, of just fitting the sigmoid
        % with new method introduced at v0.4
        if sets.edgeSettings

            spaceList = {'Otsu','Double tanh','Error function','Zscore'}; 
            
            [answer, tf] = listdlg('ListString', spaceList,...
            'SelectionMode', 'Single', 'PromptString', 'Select item', 'Initialvalue', 1,'Name', 'Make choice');

            if tf
                sets.edgeDetectionSettings.method = spaceList{answer}; 
            else
                 sets.edgeDetectionSettings.method = 'Otsu';
                import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
                sets.edgeDetectionSettings = get_default_edge_detection_settings(sets.skipDoubleTanhAdjustment);
                % user canceled or closed dialog
            end
% 
%             answer = questdlg('Which edge detection method you want to choose', ...
%             'Edge detection settings', ...
%             'Otsu','Double tanh','Error function','Zscore','Otsu');
%         
%         
%             answer = questdlg('Which edge detection method you want to choose', ...
%             'Edge detection settings', ...
%             'Otsu','Double tanh','Error function','Zscore','Otsu');
%             % Handle response
%             switch answer
%                 case 'Otsu'
%                     sets.edgeDetectionSettings.method = 'Otsu';
%                 case 'Double tanh'
%                     sets.edgeDetectionSettings.method = 'Double tanh';
%                 case 'Error function'
%                     sets.edgeDetectionSettings.method = 'Error Function';
%                 case 'Zscore'
%                     sets.edgeDetectionSettings.method = 'Zscore';
%             end

            
            
%             % add edge detection choices here
%             prompt = {'Skip double tan (default is 0, use double tan)'};
%             title = 'Edge detection settings';
%             dims = [1 35];
%             definput = {'0'};
%             answer = inputdlg(prompt,title,dims,definput);

%             sets.skipDoubleTanhAdjustment = str2double(answer{1});

            % would like to adjust these based on the experimental conditions!

        end

        if sets.random.generate
            % here we add settings
            prompt = {'Number of cut-outs','Length of a cut-out'};
            title = 'Cut-out settings';
            dims = [1 35];
            definput = {'10','50'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.random.noOfCutouts = str2double(answer{1}); % number of random cutouts from the input set
            sets.random.cutoutSize = str2double(answer{2});    
        end

        if  sets.subfragment.generate 
            prompt = {'Number of subfragments'};
            title = 'Number of subfragments settings';
            dims = [1 35];
            definput = {'2'};
            answer = inputdlg(prompt,title,dims,definput);

            sets.subfragment.numberFragments = str2double(answer{1}); % number of random cutouts from the input set
        end
        

end

