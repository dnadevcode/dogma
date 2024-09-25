function [cache] = make_consensus( lm, ts, cache )
  if nargin < 3
        cache = containers.Map();
    end

    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    flmbs3 = FancyListMgrBtnSet();
    flmbs3.NUM_BUTTON_COLS = 1;
    
    flmbs3.add_button(make_cns(ts));

    lm.add_button_sets(flmbs3);
  
    function [btnGenerateConsensus] = make_cns(ts)
        function on_make_cns(lm, ts)
            hcaSessionStruct = cache('hcaSessionStruct');
            
            
            
            % ssd align
            skipDoubleTanhAdjustment = true;
            import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
            edgeDetectionSettings = get_default_edge_detection_settings(skipDoubleTanhAdjustment);
            % 
            % import CBT.Consensus.Import.get_default_consensus_settings;
            % [barcodeConsensusSettings, dbmSettingsStruct] = get_default_consensus_settings();
            %         

            % import CBT.Hca.Import.get_barcode_params;
            % [prestretchUntrustedEdgeLenUnrounded_pixels, prestretchPixelWidth_nm,psfSigmaWidth_nm,deltaCut] = get_barcode_params(dbmSettingsStruct);

            barcodeConsensusSettings.prestretchPixelWidth_nm = 130;
            barcodeConsensusSettings.psfSigmaWidth_nm = 300;
            barcodeConsensusSettings.deltaCut = 3;
            barcodeConsensusSettings.prestretchUntrustedEdgeLenUnrounded_pixels = barcodeConsensusSettings.deltaCut * barcodeConsensusSettings.psfSigmaWidth_nm / barcodeConsensusSettings.prestretchPixelWidth_nm;

            %barcodeConsensusSettings.prestretchPixelWidth_nm = prestretchPixelWidth_nm;
            %barcodeConsensusSettings.psfSigmaWidth_nm = psfSigmaWidth_nm;
            %barcodeConsensusSettings.deltaCut = deltaCut;


            for i=1:length(hcaSessionStruct.unalignedKymos)
                kymoToAlign = double(hcaSessionStruct.unalignedKymos{i});
               tic
                [ hcaSessionStruct.alignedKymo{i},alignedKymBitMask, backgroundKym,ssdCoef] = ssd_align_first(kymoToAlign,barcodeConsensusSettings,edgeDetectionSettings);
                toc
            end


            % add option for different alignment possibilities.
%             import OptMap.KymoAlignment.NRAlign.nralign;
% 
%             for i=1:length(hcaSessionStruct.unalignedKymos)
%                 tic
%              %   hcaSessionStruct.alignedKymo{i} = nralign(hcaSessionStruct.unalignedKymos{i});
%                 toc
%             end
                      
            % this will store unfiltered barcode structure.
            hcaSessionStruct.barcodeGen=cell(length(hcaSessionStruct.unalignedKymos),1);
            
            % input pars.
            import CBT.Consensus.Import.get_default_consensus_settings;
            [barcodeConsensusSettings, dbmSettingsStruct] = get_default_consensus_settings();
%                 
            import CBT.Hca.Import.get_barcode_params;
            [prestretchUntrustedEdgeLenUnrounded_pixels, prestretchPixelWidth_nm,psfSigmaWidth_nm,deltaCut] = get_barcode_params(dbmSettingsStruct);

            barcodeConsensusSettings.prestretchUntrustedEdgeLenUnrounded_pixels = prestretchUntrustedEdgeLenUnrounded_pixels;
            barcodeConsensusSettings.prestretchPixelWidth_nm = prestretchPixelWidth_nm;
            barcodeConsensusSettings.psfSigmaWidth_nm = psfSigmaWidth_nm;
            barcodeConsensusSettings.deltaCut = deltaCut;
            
            % should we include filtered barcodes as well?
            filterSettings.filter=1;
            filterSettings.promptForfilterSettings = 1;

            % 
            import CBT.Hca.Import.get_filter_settings;
            [filterSettings] = get_filter_settings(filterSettings);
        
            import CBT.Hca.UI.Helper.gen_barcode_data;
            for i=1:length(hcaSessionStruct.unalignedKymos)
                [hcaSessionStruct.barcodeGen{i}] = gen_barcode_data(hcaSessionStruct.alignedKymo{i},barcodeConsensusSettings,filterSettings);
            end
            
            for i=1:length(hcaSessionStruct.barcodeGen)
                hcaSessionStruct.lengths(i) = length(hcaSessionStruct.barcodeGen{i}.rawBarcode);
            end
            
            rawBarcodes = cellfun(@(xx) xx.rawBarcode, hcaSessionStruct.barcodeGen,'UniformOutput',false);
            lens = cellfun(@length,rawBarcodes);

        
            if filterSettings.filter ==1
            	hcaSessionStruct.barcodeGenFiltered = cell(length(hcaSessionStruct.unalignedKymos),1);
            	filterSize = barcodeConsensusSettings.psfSigmaWidth_nm/barcodeConsensusSettings.prestretchPixelWidth_nm;

            	for i=1:length(hcaSessionStruct.unalignedKymos)
                    hcaSessionStruct.alignedKymoFiltered{i} =imgaussfilt(hcaSessionStruct.alignedKymo{i}(1:min(filterSettings.timeFramesNr,size(hcaSessionStruct.alignedKymo{i},1)),:), filterSize);     
                    [hcaSessionStruct.barcodeGenFiltered{i}] = gen_barcode_data(hcaSessionStruct.alignedKymoFiltered{i},barcodeConsensusSettings,filterSettings);
                end
                rawBarcodesFiltered = cellfun(@(xx) xx.rawBarcode, hcaSessionStruct.barcodeGenFiltered,'UniformOutput',false);

                lensFiltered = cellfun(@length,rawBarcodesFiltered);
                rawBitmasksFiltered = cell(1,length(lens));
                for i=1:length(hcaSessionStruct.barcodeGen)
                    hcaSessionStruct.lengthsFiltered(i) = length(rawBarcodesFiltered{i});
                end
            end
            
 
            
            % todo: an option to stretch barcodes to the same length for what
            % follows ..

            stretchFactors = ones(1,length(hcaSessionStruct.barcodeGen));
            untrPx = barcodeConsensusSettings.prestretchUntrustedEdgeLenUnrounded_pixels;
          
            import CBT.Bitmasking.generate_zero_edged_bitmask_row;
            rawBitmasks = cell(1,length(lens));
            for i=1:length(lens)
                rawBitmasks{i} = generate_zero_edged_bitmask_row(lens(i),round(stretchFactors(i) * untrPx));
                if filterSettings.filter == 1
                    rawBitmasksFiltered{i} = generate_zero_edged_bitmask_row(lensFiltered(i),round(stretchFactors(i) * untrPx));
                end
            end
           

            barcodeConsensusSettings.promptToConfirmTF = false;
                
            import CBT.Hca.UI.Helper.make_barcode_settings;
            [commonLength,clusterScoreThresholdNormalized,aborted] = make_barcode_settings(barcodeConsensusSettings,hcaSessionStruct.lengths);
           
           % commonLengthFiltered = ceil(mean(rawBarcodeLens));
            import CBT.Hca.UI.Helper.generate_consensus_wrapper;
            barcodeConsensusSettings.commonLength = commonLength;
            barcodeConsensusSettings.clusterScoreThresholdNormalized = clusterScoreThresholdNormalized;
            displayNames = hcaSessionStruct.names';

            if ~aborted
                [consensusStruct,rawBarcodes,rawBitmasks,displayNames,~] = generate_consensus_wrapper(rawBarcodes,rawBitmasks,barcodeConsensusSettings,hcaSessionStruct.barcodeGen,displayNames,lens,cache  )
                hcaSessionStruct.rawBarcodes  = rawBarcodes;
                hcaSessionStruct.consensusStruct = consensusStruct;
                hcaSessionStruct.rawBitmasks = rawBitmasks;
                hcaSessionStruct.displayNames = displayNames;
                if  filterSettings.filter ==1
                    displayNamesFiltered = hcaSessionStruct.names';
                    barcodeConsensusSettings.commonLength = ceil(mean(hcaSessionStruct.lengthsFiltered));
                    [consensusStructFiltered,rawBarcodesFiltered,rawBitmasksFiltered,displayNamesFiltered,~] = generate_consensus_wrapper(rawBarcodesFiltered,rawBitmasksFiltered,barcodeConsensusSettings,hcaSessionStruct.barcodeGenFiltered,displayNamesFiltered,lensFiltered,cache ,aborted );
                    hcaSessionStruct.rawBarcodesFiltered  = rawBarcodesFiltered;
                    hcaSessionStruct.consensusStructFiltered = consensusStructFiltered;
                    hcaSessionStruct.rawBitmasksFiltered = rawBitmasksFiltered;
                    hcaSessionStruct.displayNamesFiltered = displayNamesFiltered;
                end
            else
                hcaSessionStruct.consensusStruct = {};
                hcaSessionStruct.consensusStructFiltered = {};
            end
%           
% 
%             
%             if filterSettings.filter==1
%                 rawFilteredBarcodes = cellfun(@(xx) xx.rawBarcodeFiltered, hcaSessionStruct.barcodeGen,'UniformOutput',false);
%                 % new
%                  import CBT.Consensus.Core.generate_consensus_for_barcodes;
% 
%                 [consensusStructFiltered, cache] = generate_consensus_for_barcodes(rawFilteredBarcodes, displayNames,-1*ones(length(rawBarcodes),1),barcodeConsensusSettings, cache, rawBgs);
%                 hcaSessionStruct.consensusStructFiltered = consensusStructFiltered;
%                 
%                 lengths = cellfun(@length,hcaSessionStruct.consensusStructFiltered.clusterKeys);
%                 [~,b] = max(lengths);
%                 key = hcaSessionStruct.consensusStructFiltered.clusterKeys{b};
%                 consSt = hcaSessionStruct.consensusStructFiltered.barcodeStructsMap(key);
%                 hcaSessionStruct.names{length(lens)+1} = key ;
%                 rawFilteredBarcodes{length(lens)+1} = consSt.barcode;
%                 rawFilteredBitmasks{length(lens)+1} = logical(consSt.indexWeights);
%                 %hcaSessionStruct.consensusFiltered =  consensusStructFiltered.barcodeStructsMap(consensusStructFiltered.finalConsensusKey);
% 
%                 hcaSessionStruct.rawFilteredBarcodes = rawFilteredBarcodes;
% 
%             end
            hcaSessionStruct.filterSettings = filterSettings;
            hcaSessionStruct.barcodeConsensusSettings = barcodeConsensusSettings;
            cache('hcaSessionStruct') = hcaSessionStruct;            
        end

        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnGenerateConsensus = FancyListMgrBtn(...
            'Make consensus structure', ...
            @(~, ~, lm) on_make_cns(lm, ts));
    end
end

