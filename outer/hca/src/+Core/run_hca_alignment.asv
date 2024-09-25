function [barcodeGenC, consensusStruct, comparisonStruct, theoryStruct, hcaSets] = run_hca_alignment(hcaSets)
    %   run_hca_theory - calculates HCA theory

    %   Args: hcaSets - hca settings, either passed through
    %   hca_barcode_alignment or loaded as paramenter                

    %   Returns:
    %       barcodeGenC - barcode structure
    %       comparisonStruct - comparison structure
    %       theoryStruct - theory structure

    %   Exmple:
    %   [hcaSets,hcaaligner.names] = Core.Default.read_default_sets('hcaalignmentsets.txt');
    %   hcaSets.default.kymosets.filenames = 'kymofolder';
    %  [barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, hcaSets] = run_hca_alignment(hcaaligner.sets.default)




    % change names to support previous hca
    hcaSets.filterSettings.filter = hcaSets.filterSettingsFilter;
    hcaSets.random.generate = hcaSets.randomGenerate;
    hcaSets.subfragment.generate = hcaSets.subfragmentGenerate;

    import CBT.Hca.Settings.non_default_settings;
    hcaSets = non_default_settings(hcaSets);


%     [,hcaSets.kymosets.filenames,hcaSets.kymosets.fileext] = cellfun(@(x) fileparts(x),hcaSets.kymofolder,'un',false);
    hcaSets.kymosets.filenames = hcaSets.kymofolder;
    hcaSets.kymosets.kymofilefold = cell(1,length(   hcaSets.kymosets.filenames));
    hcaSets.output.matDirpath = fileparts(hcaSets.kymofolder{1});
    % timestamp for the results
    hcaSets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

    %% Add a choice to load barcodeGen directly
    [~,~,en] = fileparts(hcaSets.kymofolder{1});
    switch en
        case '.tif'
            % add kymographs
            import CBT.Hca.Import.add_kymographs_fun;
            [kymoStructs] = add_kymographs_fun(hcaSets);
        
            %  put the kymographs into the structure
            import CBT.Hca.Core.edit_kymographs_fun;
            kymoStructs = edit_kymographs_fun(kymoStructs,hcaSets.timeFramesNr);
        
            % align kymos
            import CBT.Hca.Core.align_kymos;
            [kymoStructs] = align_kymos(hcaSets,kymoStructs);
            
            % generate barcodes
            import CBT.Hca.Core.gen_barcodes;
            barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, hcaSets);
        case '.mat'
            bgSessionData = load(hcaSets.kymofolder{1});
            if isfield(bgSessionData,'barGenMerged')
                barcodeGen = bgSessionData.barGenMerged;
            else
                if isfield(bgSessionData,'theoryGen')
                    barcodeGen = cell(1,length(bgSessionData.theoryGen.theoryBarcodes));
                    for i=1:length(bgSessionData.theoryGen.theoryBarcodes)
                        barcodeGen{i}.rawBarcode = bgSessionData.theoryGen.theoryBarcodes{i};
                        barcodeGen{i}.rawBitmask = ones(1,length( barcodeGen{i}.rawBarcode));
                    end
                else
                barcodeGen = bgSessionData.barcodeGen;
                end
            end
            if ~isfield(barcodeGen{1}, 'name')
                for i=1:length(barcodeGen)
                    barcodeGen{i}.name = [num2str(i), 'renamed_bar'];
                end
            end
        otherwise
    end

    % shrink find
%     import Core.shrink_finder_fun;
%     shrink_finder_fun(kymoStructs, hcaSets)


    
    % generate consensus
    import CBT.Hca.Core.gen_consensus;
    consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGen,hcaSets);
    %
    % Alternative:
    %         import CBT.Hca.UI.Helper.select_all_consensus;
    %         [consensusStruct,sets] = select_all_consensus(consensusStructs,sets);
    
    % select consensus
    import CBT.Hca.UI.Helper.select_consensus
    [consensusStruct,hcaSets] = select_consensus(consensusStructs,hcaSets);


  % generate random cut-outs
    if hcaSets.random.generate
        % % Here generate cutouts using M_files_HCA_exp_cutouts:
        import CBT.Hca.Core.Random.cutout_barcodes;
        barcodeGenRandom = cutout_barcodes(barcodeGen,hcaSets);
        barcodeGenC = barcodeGenRandom;
    else
        barcodeGenC = barcodeGen;
    end

    if  hcaSets.subfragment.generate
        import CBT.Hca.Core.Subfragment.gen_subfragments;
        barcodeGenC = gen_subfragments(barcodeGen,hcaSets.subfragment.numberFragments); 
    else
        barcodeGenC = barcodeGen;
    end


    if  isempty(hcaSets.theoryfolder{1})
        % load(thryFile);
        [hcaSets.theoryFile{1},hcaSets.theoryFileFold{1}] = uigetfile(pwd,strcat(['Select single theory .mat file to process']));
    else
        for i=1:length(hcaSets.theoryfolder)
            [hcaSets.theoryFileFold{i},mid,en] = fileparts(hcaSets.theoryfolder{i});
            hcaSets.theoryFile{i} = [mid,en];
        end
    end
    hcaSets.theory.precision = 5;
    hcaSets.theory.theoryDontSaveTxts = 1;
    import CBT.Hca.UI.Helper.load_theory;
    theoryStruct = load_theory(hcaSets);

    
    % extract from name
    hcaSets.theory.nmbp = hcaSets.nmbp;
    hcaSets.theory.stretchFactors = 1-hcaSets.sF/100:hcaSets.sFstep/100:1+hcaSets.sF/100; % length re-scaling factor

    import CBT.Hca.Core.Analysis.convert_nm_ratio;
    theoryStruct = convert_nm_ratio(hcaSets.theory.nmbp, theoryStruct,hcaSets );
    
    
    hcaSets.theoryFile=[];
    hcaSets.theoryFileFold = [];
    
    % compare theory to experiment
    import CBT.Hca.Core.Comparison.compare_distance;
    [rezMax,bestBarStretch,bestLength] = compare_distance(barcodeGenC,theoryStruct, hcaSets, consensusStruct );

    comparisonStructAll = rezMax;
    for i=1:length(comparisonStructAll)
        for j=1:length(bestBarStretch{i})
            comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
            comparisonStructAll{i}{j}.length = bestLength{i}(j);
        end
    end
    import CBT.Hca.Core.Comparison.combine_theory_results;
    [comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);
    % 
    % 
    % 
    %    % assign all to base
    if hcaSets.saveOutput
            assignin('base','barcodeGenC', barcodeGenC)
            assignin('base','theoryStruct', theoryStruct)
            assignin('base','hcaSets', hcaSets)
            assignin('base','rezMax', rezMax)
    end
        
%             import CBT.Hca.UI.Helper.plot_any_bar;
%             plot_any_bar(1,barcodeGenC,rezMax,theoryStruct,1);

        import CBT.Hca.Core.additional_computations
        additional_computations(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct,comparisonStructAll,hcaSets );

        if hcaSets.identifyDisc
            import Core.identify_discriminative;
            [is_distinct,numMatchingSpecies,uniqueMatchSequences,refNums,refNumBad] = identify_discriminative(hcaSets,barcodeGenC,rezMax, {theoryStruct.name});
            
            for i=1:length(comparisonStruct)
                comparisonStruct{i}.discriminative.is_distinct = is_distinct(i);
                comparisonStruct{i}.discriminative.numMatchingSpecies = numMatchingSpecies(i);
                comparisonStruct{i}.discriminative.uniqueMatchSequences = uniqueMatchSequences(i);
            end
            if hcaSets.saveOutput
                    assignin('base','is_distinct', is_distinct);
            end
        end


%               import CBT.Hca.UI.get_display_results;
%         get_display_results(barcodeGenC,consensusStruct, comparisonStruct, theoryStruct, sets);

end

