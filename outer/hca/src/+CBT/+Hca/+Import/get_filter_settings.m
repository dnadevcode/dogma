function [filterSettings] = get_filter_settings(filterSettings)
    %  get_filter_settings
    
    titleText = 'Selection of number of timeframes for filtered kymographs';

    if (filterSettings.promptForfilterSettings == 0)
        filterSettings.filter = 1;
        filterSettings.timeFramesNr =1;
        return 
    end

    % Set a threshold for grouping barcodes.
    continueGenChoice = questdlg('Add filtered barcodes', 'Add filtered barcodes?', 'Yes (recommended)', 'No', 'Yes (recommended)');
    filterSettings.filter = strcmp(continueGenChoice, 'Yes (recommended)');

    if filterSettings.filter
        % this chooses timeframe nr for filtered

        % choose timeframes nr for unfiltered kymographs
        import CBT.Hca.UI.get_hca_settings;
        [ filterSettings.timeFramesNr] = get_hca_settings(filterSettings.timeFramesNr,titleText); 

        % choose timeframes nr for unfiltered kymographs
        import CBT.Hca.UI.get_hca_filter_settings;
        [ filterSettings.filterSize] = get_hca_filter_settings(filterSettings.filterSize,'select filter size'); 

        % choose when to filter          
        filterChoice = questdlg('Filter before generating barcodes ', 'Choose when to filter', 'before (recommended)', 'after', 'before (recommended)');
        filterSettings.filterMethod = strcmp(filterChoice, 'before (recommended)');
    end

end