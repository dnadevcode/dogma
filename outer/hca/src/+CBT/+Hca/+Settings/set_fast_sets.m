function [ sets ] = set_fast_sets(adddefaulttheorysets)
    % set_fast_sets 
    % Set quick pre-defined settings for HCA.
    % If we want, we can run the program
    % without any prompts, good if we want to run multiple things.
    %
    %     Returns:
    %         sets (struct): pre-defined settings structure
    % 
    %     Example:
    %         This is an example: run sets = set_fast_sets()
    
    if nargin < 1
        adddefaulttheorysets = 0;
    end
    %
    sets.kymosets.askforkymos = 1; % sould we ask for kymos
    
    % sample filename, comment out for the release version
    sets.kymosets.readfirstfolder = 0; % should we read all the tifs from current folder for this
   % sets.kymosets.kymofilefold{1} = '/home/albyback/rawData/dnaData/humanData/exp/2330P18/Raw Kymos/';
   sets.kymosets.kymofilefold{1} = '';
   sets.kymosets.filenames{1} = [];
    
    % all the different setting choices
    sets.kymosets.askforsets = 1; % should we ask for settings
	sets.timeFramesNr = 0; % 0 - take all timeframes
    sets.alignMethod = 1; % 1 - nralign, 2 - ssdalign, 3 - ssg
   	sets.filterSettings.filter = 0; % filter
    sets.genConsensus = 0; % generate consensus
    sets.random.generate = 0; % generate random

    
    % filter settings. For one timeframe, refer to P.Torche paper for best
    % parameters
	sets.filterSettings.filterMethod = 0; % 0 - filter after stretching, 1 - before
    sets.filterSettings.filterSize = 2.3;
    
    % edge detection  settings
    sets.skipEdgeDetection = 0; % skip edge detection
    sets.skipDoubleTanhAdjustment = 1;
    import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
    sets.edgeDetectionSettings = get_default_edge_detection_settings(sets.skipDoubleTanhAdjustment);
    
    % bitmask settings
    sets.bitmasking.prestretchPixelWidth_nm = 130; % camera nm/px
	sets.bitmasking.psfSigmaWidth_nm = 300; % psf in nanometers
	sets.bitmasking.deltaCut = 3; % how many delta's should we take
    sets.bitmasking.untrustedPx = sets.bitmasking.deltaCut*sets.bitmasking.psfSigmaWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;

    % consensus settings
    sets.consensus.barcodeNormalization = 'bgmean';
    sets.consensus.promptForBarcodeClusterLimit = 1;
    sets.consensus.threshold = 0.75;

    % add comparison with DTW:
    sets.comparison.useDTW = 1;
    % add comparison with PCC:
    sets.comparison.usePCC = 1;

    
    % Random cutouts
    sets.random.noOfCutouts = 2; % number of random cutouts from the input set
    sets.random.cutoutSize = 200; % size of region to be cut out (units: pixels)

    %
    sets.theory.askfortheory = 1;
% 	sets.theoryFileFold{1} = '/home/albyback/rawData/dnaData/humanData/';
%     sets.theoryFile{1}= 'ncbi_chromosomes_181022_130nmPixel_0.3nmPERbp.mat';
    sets.theoryFileFold{1} = '';
    sets.theoryFile{1} = '';
    sets.theory.askfornmbp = 1;
    sets.theory.nmbp = 0.225;
	%sets.theory.skipStretch = 1; %0 - do not stretch,  1 - stretch
	sets.theory.stretchFactors = [0.95 0.96 0.9700    0.9800    0.9900    1.0000    1.0100    1.0200    1.0300 1.04 1.05];
%         % 
    
    sets.export.savetxt = 1;
    sets.displayResults = 1;
    
    
   % sets.pvalue.fold = '/home/albyback/git/Projects/hca/PVAL/';
    sets.pvalue.file = ''; 
	sets.pvalue.promtforparams = 0;
    sets.pvalue.len2 = 5000000;
    sets.pvalue.pixelWidth_nm = 130;
    sets.pvalue.lenMin = 100;
    sets.pvalue.lenMax = 1000;
    sets.pvalue.numRnd = 1000;
    sets.pvalue.psfSigmaWidth_nm = 300;
    
  %  sets.output.matDirpath =  '/home/albyback/rawData/dnaData/humanData/output/';
    sets.output.matDirpath = '';
    sets.output.askforoutputdir = 1;
    % sets.export.savetxt = 1;
%     sets.displayResults = 1;
%     
% 
% 
% 
%  
%     % stretching
%     sets.skipPrechoice = 1;
%     sets.prestretchMethod = 0; % 0 - do not prestretch % 1 - prestretch to common length

%     % filter settings
%     sets.skipFilterSettings = 1;
%     
%    
%     % Theory generation
%     sets.skipBarcodeGenSettings = 1;
%     sets.skipChangeBpNmRatio = 1;
% 
%     % all file folders should be possibly removed in the release

% 
%     % barcode theory generation settings
%     sets.barcodeGenSettings.meanBpExt_nm = 0.225;
%     sets.barcodeGenSettings.pixelWidth_nm = 130;
%     sets.barcodeGenSettings.concNetropsin_molar = 6;
%     sets.barcodeGenSettings.concYOYO1_molar = 0.02;
%     sets.barcodeGenSettings.isLinearTF = 0;
%     sets.barcodeGenSettings.psfSigmaWidth_nm = 300;
%     sets.barcodeGenSettings.deltaCut = 3;
%     sets.barcodeGenSettings.widthSigmasFromMean = 4;
%     sets.barcodeGenSettings.yoyo1BindingConstant = 26;
%     
%     sets.skipNullModelChoice = 1;
%     sets.nullModelPath = '/home/albyback/git/WORKSHOP/HCA_v1.6/nullmodel';
%     sets.askForPvalueSettings = 0;
%     sets.pvaluethresh = 0.01;
%     sets.contigSettings.numRandBarcodes = 1000;

    if adddefaulttheorysets
        sets.theoryGen.meanBpExt_nm = 0.3;
        sets.theoryGen.pixelWidth_nm = 130;
        sets.theoryGen.concN = 6;
        sets.theoryGen.concY = 0.02;
        sets.theoryGen.concDNA = 0.2;
        sets.theoryGen.psfSigmaWidth_nm = 300;
        sets.theoryGen.deltaCut = 3;
        sets.theoryGen.isLinearTF = 0;
        sets.theoryGen.widthSigmasFromMean = 4;
        sets.theoryGen.computeFreeConcentrations = 1;
        sets.theoryGen.model = 'literature';
        sets.lambda.fold = 'tests/data/';
        sets.lambda.name = 'sequence.fasta';
        import CBT.Hca.Core.Theory.cb_model;
        sets.model = cb_model();
    end

end

