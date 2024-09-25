function [theoryGen] = run_hca_theory(hcaSets,yoyoConst, netropsinConst,customModel)
    %   run_hca_theory - calculates HCA theory

    %   Args: hcaSets
    %           yoyoConst - in case user want to vary the YOYO-1 binding constant
    %           netropsinConst - in case one wants to 
    if nargin <2
        yoyoConst = 26;
        netropsinConst = 0.4;
    end


    % make sets compatible with prev structure
    hcaSets.theoryGen.method = hcaSets.method;
    hcaSets.theoryGen.computeFreeConcentrations =  hcaSets.computeFreeConcentrations;
    hcaSets.theoryGen.concN = hcaSets.concN;
    hcaSets.theoryGen.concY = hcaSets.concY;
    hcaSets.theoryGen.concDNA = hcaSets.concDNA;
    hcaSets.theoryGen.psfSigmaWidth_nm = hcaSets.psfSigmaWidthNm;
    hcaSets.theoryGen.pixelWidth_nm = hcaSets.pixelWidthNm;
    hcaSets.theoryGen.meanBpExt_nm = hcaSets.meanBpExtNm;
    hcaSets.theoryGen.isLinearTF = hcaSets.isLinearTF;

    hcaSets.lambda.fold  =  hcaSets.fold ;
    hcaSets.lambda.name  =  hcaSets.name ;
    hcaSets.theoryGen.k =  max(2.^15,2.^(hcaSets.k));
    hcaSets.theoryGen.m  = min(2.^15,2.^(hcaSets.m));
    hcaSets.theoryGen.computeBitmask = hcaSets.computeBitmask;
    % timestamp to add to theories name
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');


    import CBT.Hca.Core.Theory.choose_cb_model;
    [hcaSets.model ] = choose_cb_model(hcaSets.theoryGen.method,hcaSets.pattern, yoyoConst, netropsinConst);

    if nargin >=4
        hcaSets.model.netropsinBindingConstant = customModel;
    end

    % theories names
    theories = hcaSets.folder;

    try
        hcaSets.resultsDir = fullfile(fileparts(theories{1}),'theoryOutput'); % we always save theories output in the same folder as provided data
        % make theoryData folder
        [~,~] = mkdir(hcaSets.resultsDir);

    catch
        hcaSets.resultsDir ='';
    end
    

    % compute free concentrations
    import CBT.Hca.Core.Theory.compute_free_conc;
    hcaSets = compute_free_conc(hcaSets);
    
    theoryGen = struct();

    % save sets
    theoryGen.sets = hcaSets.theoryGen;
    meanBpExt_nm = hcaSets.meanBpExtNm;
    pixelWidth_nm = hcaSets.pixelWidthNm;
    psfSigmaWidth_nm = hcaSets.psfSigmaWidthNm;
    linear = hcaSets.isLinearTF;

    theoryBarcodes = cell(1,length(theories));
    theoryBitmasks = cell(1,length(theories));
    theoryNames = cell(1,length(theories));
    theoryIdx = cell(1,length(theories));
    bpNm = cell(1,length(theories));

    import CBT.Hca.Core.Theory.compute_theory_barcode;

    if length(theories) > 1

    % loop over theory file folder
    parfor idx = 1:length(theories)
        try
            addpath(genpath(fileparts(theories{idx})));
            disp(strcat(['loaded theory sequence ' theories{idx}] ));
        catch
            addpath(genpath(fileparts(theories{idx}.Filename))); % if preloaded

        end
    
        % new way to generate theory, check theory_test.m to check how it works
        [theorySeq, header,bitmask] = compute_theory_barcode(theories{idx},hcaSets);
    
	    theoryBarcodes{idx} = theorySeq;
        theoryBitmasks{idx} = bitmask;
    
        theoryNames{idx} = header;
        theoryIdx{idx} = idx;
        bpNm{idx} = meanBpExt_nm/psfSigmaWidth_nm;         
    end
    else
    % loop over theory file folder
    for idx = 1:length(theories)
        try
            addpath(genpath(fileparts(theories{idx})));
            disp(strcat(['loaded theory sequence ' theories{idx}] ));
        catch
            addpath(genpath(fileparts(theories{idx}.Filename))); % if preloaded
        end
    
        % new way to generate theory, check theory_test.m to check how it works
        [theorySeq, header,bitmask] = compute_theory_barcode(theories{idx},hcaSets);
    
	    theoryBarcodes{idx} = theorySeq;
        theoryBitmasks{idx} = bitmask;
    
        theoryNames{idx} = header;
        theoryIdx{idx} = idx;
        bpNm{idx} = meanBpExt_nm/psfSigmaWidth_nm;         
    end

    end
    theoryGen.theoryBarcodes = theoryBarcodes;
    theoryGen.theoryBitmasks = theoryBitmasks;
    theoryGen.theoryNames = theoryNames;
    theoryGen.theoryIdx = theoryIdx;
    theoryGen.bpNm = bpNm;
    
    matFilename = strcat(['theoryGen_', num2str(meanBpExt_nm) '_' num2str(pixelWidth_nm) '_' num2str(psfSigmaWidth_nm) '_' num2str(linear) '_' sprintf('%s_%s', timestamp) 'session.mat']);
    matFilepath = fullfile(hcaSets.resultsDir, matFilename);
    
    disp('Finished calculating theories')
    if nargout < 1
        assignin('base','theoryGen', theoryGen)
        disp('Assigned theoryGen to workspace');
        
        
        save(matFilepath, 'theoryGen','-v7.3');
        fprintf('Saved theory mat filename ''%s'' to ''%s''\n', matFilename, matFilepath);
    end

end

