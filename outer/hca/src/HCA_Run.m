function [hcaStruct] = HCA_Run(sets, hcaStruct)
    % HCA_Run
    % Used for comparing fagments of human chromosome to chromosomes using CB (competitive binding) theory   
    %
    %     Args:
    %         sets (struct): Input settings to the method
    %         hcaStruct (struct): Input structure, if non-empty, load
    %         result structure instead of computing everything
    % 
    %     Returns:
    %         hcaStruct: Return structure
    % 
    %     Example:
    %         This is an example: run [hcaStruct] = HCA_Gui(sets) 
    %
    
    % TODO: return warnings to the places where it is important to know
    % what is being done and to check for consistency
    
    % timestamp for the results
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

    if nargin < 1 % if settings were not provided
        % import settings
        import CBT.Hca.Import.import_hca_settings;
        [sets] = import_hca_settings('hca_settings.txt');
    end
    
    if nargin < 2 
        % load data options: use load_data wrapper after selecting
        % settings.

        % first option: regular kymographs
        
        % load all user selected settings
        import CBT.Hca.Settings.get_user_settings;
        sets = get_user_settings(sets);
        sets.timestamp = timestamp;
        mkdir(sets.timestamp);

        % data wrapper after user settings
        try
            % check if it's a mat file (consensus/theory)
            % or kymographs (.tif)
            [ff,fd,fl] = fileparts(sets.kymosets.filenames{1});
        catch
            fl = 'rand';
        end
        
        switch fl
            case '.mat'
                disp("Running consensus vs theory comparison, some options i.e. subfragment generation are skipped");
                sets.type = 'consensus';
                import CBT.Hca.Import.load_data_wrapper;
                [uniqueNames,uniqueFolds,barcodeGenC] = load_data_wrapper(sets.type,sets.output.matDirpath,sets,sets.kymosets.kymofilefold{1} );
            case '.tif'
                % add kymographs
                import CBT.Hca.Import.add_kymographs_fun;
                [kymoStructs] = add_kymographs_fun(sets);

                %  put the kymographs into the structure
                import CBT.Hca.Core.edit_kymographs_fun;
                kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr);

                % align kymos
                import CBT.Hca.Core.align_kymos;
                [kymoStructs] = align_kymos(sets,kymoStructs);

                % generate barcodes
                import CBT.Hca.Core.gen_barcodes;
                barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);

                % generate consensus
                import CBT.Hca.Core.gen_consensus;
                consensusStructs = CBT.Hca.Core.gen_consensus(barcodeGen,sets);
                %
                    
                % change the select_consensus ... And then save in a simple
                % way
                % select consensus
                
%                 select_same_cut_consensus
                import CBT.Hca.UI.Helper.select_same_cut_consensus
                [consensusStruct,sets] = select_same_cut_consensus(consensusStructs,sets);
% 
%                 figure,hist(cellfun(@(x) length(x.rawBarcode),barcodeGen(consensusStruct{1}.barcodesInConsensus)))
% 
%                 figure,histogram(cellfun(@(x) length(x.rawBarcode),barcodeGen),'Normalization','count','BinMethod','integers')
%                 hold on
%                 
%                 lengths = cellfun(@(y) cellfun(@(x) length(x.rawBarcode),barcodeGen(y.barcodesInConsensus)),consensusStruct,'un',false);
%                 histogram(lengths(1),'Normalization','count','BinMethod','integers')
%                 histogram(lengths(2),'BinMethod','integers')
%                 histogram(lengths(3),'BinMethod','integers')
%                 legend({'All lengths histogram','Clust1','Clust2','Clust3'})

%                 import CBT.Hca.UI.Helper.select_consensus
%                 [consensusStruct,sets] = select_consensus(consensusStructs,sets);

                    %% These two options I am not sure if I should keep here?
                % generate random cut-outs
                if sets.random.generate
                    % % Here generate cutouts using M_files_HCA_exp_cutouts:
                    import CBT.Hca.Core.Random.cutout_barcodes;
                    barcodeGenRandom = cutout_barcodes(barcodeGen,sets);
                    barcodeGenC = barcodeGenRandom;
                else
                    barcodeGenC = barcodeGen;
                end

                if  sets.subfragment.generate
                    import CBT.Hca.Core.Subfragment.gen_subfragments;
                    barcodeGenC = gen_subfragments(barcodeGen,sets.subfragment.numberFragments); 
                else
                    barcodeGenC = barcodeGen;
                end

            otherwise
                    % if nothing selected, just create synthetic data, with the
                % following parameters 
                disp("generating random data since no familiar data selected");
                sets.type = 'synthetic';
                sets.kymosets.kymofilefold{1} = [];
                sets.svType = 1; % type of sv, insertion
                sets.numSamples = 10; % number of random samples
                sets.length1 =500; % len1 
                sets.length2 = 50; % len2
                sets.circ = 1; % whether is circular
                sets.kernelSigma = 2.3; % kernel sigma
                sets.pccScore = 0.9; % pcc score
                sets.c = sets.w; % c, if mp is used
                sets.genConsensus = 0;
                import CBT.Hca.Import.load_data_wrapper;
                [uniqueNames,uniqueFolds,barcodeGenC] = load_data_wrapper(sets.type,sets.output.matDirpath,sets,sets.kymosets.kymofilefold{1} );

        end
         

        % all these things are processing of input data. But could still
        % be in load data wrapper as they are not really producing any
        % results just yet
        %         
        % TODO : add consensus plot
        %         % recreate consensus result from cbc_gui? 
        %         import CBT.Hca.Core.gen_consensus_plot;
        %         gen_consensus_plot(consensusStructs,sets);

        % this is really processing data rather than loading so should be here


    
        %% now user theories. They could already be in txt files (if generated 
        % (with HCA 4.0.0), but we should keep support for older theory files too 
        % Also, this could be generated via load data wrapper!
        sets.theoryFile=[];
        sets.theoryFileFold = [];

        % get user theory
        import CBT.Hca.Settings.get_user_theory;
        [theoryStruct, ~] = get_user_theory(sets);

        % split the comparison into regular & consensus (of course
        % depending on the type of data that we have, if it is rawbarcodes
        % + consensus, or just rawbarcodes, or just consensus.
        % in the barcodeGen, rawBarcode, rawBitmask, consensusBarcode,
        % consensusBitmask should be the structure for these, then it is
        % easy to differentiate  if we're using "raw" or "consensus", also
        % if we're running "singleframe", when i.e. one kymograph is split
        % into all of it's single frames.
        sets.genConsensus = 0;
        % compare theory to experiment
         import CBT.Hca.Core.Comparison.compare_distance;
         [rezMax,bestBarStretch,bestLength] = compare_distance(barcodeGenC,theoryStruct, sets, [] );
%         import CBT.Hca.Core.Comparison.compare_theory_to_exp;
%         comparisonStructAll = compare_theory_to_exp(barcodeGenC, theoryStruct, sets, consensusStruct);
        
    % change notation a little from previous, shouldn't change the speed
    % too much
      comparisonStructAll = rezMax;
        for i=1:length(comparisonStructAll)
            for j=1:length(bestBarStretch{i})
                comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
                comparisonStructAll{i}{j}.length = bestLength{i}(j);
            end
        end
        import CBT.Hca.Core.Comparison.combine_theory_results;
        [comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);
        
        if max(cellfun(@(x) x.maxcoef(1),comparisonStruct)) == 0
            disp("No barcodes had length longer than preselected window length");
            return;
        end

%         % compare theory to experiment
%         import CBT.Hca.Core.Comparison.compare_theory_to_exp;
%         comparisonStructAll = compare_theory_to_exp(barcodeGenC, theoryStruct, sets, consensusStruct);


%         % bugcheck: if only one theory
%         import CBT.Hca.UI.combine_chromosome_results;
%         comparisonStruct = combine_chromosome_results(theoryStruct,comparisonStructAll);

        sets.userDefinedSeqCushion = 50; % move to sets
       % assign all to base
        assignin('base','barcodeGenC', barcodeGenC)
        assignin('base','theoryStruct', theoryStruct)
        assignin('base','comparisonStructAll', comparisonStructAll)
        assignin('base','sets', sets)
        assignin('base','comparisonStruct', comparisonStruct)
        assignin('base','consensusStruct', [])

        %
      %  sets.displayResults = 1;
        sets.plotallmatches = 0;
        import CBT.Hca.UI.get_display_results;
        get_display_results(barcodeGenC,[], comparisonStruct, theoryStruct, sets);


        import CBT.Hca.Core.additional_computations
        additional_computations(barcodeGenC,[], comparisonStruct, theoryStruct,comparisonStructAll,sets );
        
        disp("Dist scores calculates and results plotted and saved, now checking for p-values...");
        % whether to compute p-values ? 
%         sets.computepval = 1;
        
        import CBT.Hca.Core.Pvalue.compute_pvals_simple;
        [pval] = compute_pvals_simple(theoryStruct,comparisonStruct,sets);

        % If one computes these directly from data, then we want a random
        % theory and 1000 random barcodes. We then compare compare these (
        % so takes 1000 times) using compare distance. Finally we compute
        % p-values
     
      
      %  sets.output.matDirpath = '/home/albyback/rawData/dnaData/humanData/output/';
  
    else
        % run load results function or something similar..
    end
    
end
