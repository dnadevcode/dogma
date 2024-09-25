function [cache] = export_consensus_structure( lm, cache )
    % make_consensus_structure
    % input lm, cache
    if nargin < 2
        cache = containers.Map();
    end

	% creates a button set
    import Fancy.UI.Templates.create_button_set;
    create_button_set(lm, @make_cns);
    
    function [btnGenerateConsensus] = make_cns()
        function on_make_cns()
            % load structures from cache
            hcaSessionStruct = cache('hcaSessionStruct');
            sets = cache('sets');  
            hcaSessionStruct.sets = sets;
            
            % export usual structure
            import CBT.Hca.Export.export_results_mat;
            export_results_mat(hcaSessionStruct, 'barcode_struct');
            
            % also, if there is consensus, maybe we want to export in a
            % structure that is loadable with CBC:
            import CBT.Hca.Export.export_cbc_compatible_consensus;
            export_cbc_compatible_consensus(hcaSessionStruct);

        end

        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnGenerateConsensus = FancyListMgrBtn(...
            'Make consensus structure', ...
            @(~, ~, lm) on_make_cns());
    end
end

