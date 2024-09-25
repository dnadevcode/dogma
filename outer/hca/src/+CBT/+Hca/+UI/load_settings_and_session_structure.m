function [cache] = load_settings_and_session_structure( lm, cache )
    % load_settings_and_session_structure
    % Loads settings and session structure
    %
    %     Args:
    %         lm (list): List manager
    %         cache (container): cached results
    % 
    %     Returns:
    %         cache (container): cached results
    %

    if nargin < 2
        cache = containers.Map();
    end

    % creates a button set
    import Fancy.UI.Templates.create_button_set;
    create_button_set(lm, @edit_kymographs);
    
    % edit_kymographs sets up the structure
    function [btnGenerateConsensus] = edit_kymographs()
        function on_edit_kymographs(lm) 
            % load default settings. Maybe let user select default settings
            % file?
            disp('Loading default settings file, CBT.Hca.Import.set_default_settings');
            import CBT.Hca.Import.set_default_settings;
            sets = set_default_settings();
            
            % choose timeframes nr for unfiltered kymographs. In case
            % number of timeframes is set to 0, we take all the available
            % timeframes
            import Fancy.UI.Templates.get_frame_nr;
            [ sets.timeFramesNr] = get_frame_nr(sets.promptForTimeFr,sets.timeFramesNr,'Selection of number of timeframes'); 
            
            % we extract kymo's from a list
            import Fancy.UI.Templates.extract_kymos_from_list;
            [hcaSessionStruct] = extract_kymos_from_list(lm, sets.timeFramesNr);
            
            % put in the cache for future reference
            cache('hcaSessionStruct') = hcaSessionStruct;     
            cache('sets') = sets;     
        end
        
        buttonName = 'Load settings and add kymographs to struture';
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnGenerateConsensus = FancyListMgrBtn(buttonName, @(~,~, lm) on_edit_kymographs(lm));
    end


end



