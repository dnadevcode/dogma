function [cache] = make_consensus_structure( lm, cache )
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
            
            % choose alignment method to align the kymographs
            if ~sets.skipAlignChoice
                alignChoice = questdlg('Choose alignment method', 'Two possible alignment methods', 'nralign', 'ssdalign',  'ssdalign');
                sets.alignMethod = ~strcmp(alignChoice,'nralign');
            end
            
            % kymo alignment. We run this to align kymographs
            import CBT.Hca.UI.align_kymos;
            [hcaSessionStruct.alignedKymo,hcaSessionStruct.leftEdgeIdxs,hcaSessionStruct.rightEdgeIdxs] = align_kymos(sets,hcaSessionStruct.unalignedKymos);
            
            % this will store unfiltered barcode structure.
            hcaSessionStruct.barcodeGen = cell(length(hcaSessionStruct.unalignedKymos),1);
            
            % input parameters.
            import CBT.Hca.Import.get_barcode_params;
            sets.barcodeConsensusSettings = get_barcode_params(sets.barcodeConsensusSettings,sets.skipDefaultConsensusSettings);
            
            % here we make a choice for filtering the barcodes
            import CBT.Hca.Import.get_filter_settings;
            [sets.filterSettings] = get_filter_settings(sets.filterSettings);
        
            % choose pre-stretching method
            if ~sets.skipPrechoice
                prestretchChoice = questdlg('Prestretch barcodes to the same lengths', 'Prestretching', 'yes (if multiple copies of one DNA molecule)', 'no', 'yes (if multiple copies of one DNA molecule)');
                sets.prestretchMethod = strcmp(prestretchChoice, 'yes (if multiple copies of one DNA molecule)');
            end
            
            % generate barcodes
            import CBT.Hca.UI.Helper.gen_barcodes;
            hcaSessionStruct = gen_barcodes(hcaSessionStruct,sets);
            
            % what does this do?
            sets.barcodeConsensusSettings.promptToConfirmTF = false;
                
            if ~sets.skipbarcodeConsensusSettings
                import CBT.Hca.UI.Helper.make_barcode_settings;
                [commonLength,aborted] = make_barcode_settings(hcaSessionStruct.lengths);
                sets.barcodeConsensusSettings.aborted = aborted;
                sets.barcodeConsensusSettings.commonLength = commonLength;
               % sets.barcodeConsensusSettings.barcodeClusterLimit = clusterScoreThresholdNormalized;
            end

            % generate consensus
            import CBT.Hca.UI.Helper.gen_consensus
            hcaSessionStruct = gen_consensus(hcaSessionStruct,sets);
            
            % select consensus
            import CBT.Hca.UI.Helper.select_consensus
            [hcaSessionStruct,sets] = select_consensus(hcaSessionStruct,sets);

            cache('hcaSessionStruct') = hcaSessionStruct;     
            cache('sets') = sets;            
        end

        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnGenerateConsensus = FancyListMgrBtn(...
            'Make consensus structure', ...
            @(~, ~, lm) on_make_cns());
    end
end

