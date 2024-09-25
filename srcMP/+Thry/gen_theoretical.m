function [theoryStruct,theoryGen,barcodeGen] = gen_theoretical(fastaFile,meanBpExt_nm,isLinearTF,pixelWidth_nm,setFile)


% fastaFile =  '/home/albyback/git/bargroupingprototype/test/DA32087.fasta';



% timestamp for the results
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

% import settings
if nargin < 5
    setFile = 'theory_settings.txt';
end
import CBT.Hca.Import.import_settings;
[sets] = import_settings(setFile);


sets.fastas = 'fastas.txt';
fd = fopen(sets.fastas,'w');
for i=1:length(fastaFile)
    fprintf(fd,'%s\n',fastaFile{i});
end
fclose(fd);

sets.skipBarcodeGenSettings = 1;
sets.promptfortheory = 0;
sets.promptforsavetheory = 0;
sets.savetxts = 0;
% load default settings
import CBT.Hca.Settings.get_user_theory_sets;
sets = get_user_theory_sets(sets);


if ~sets.skipBarcodeGenSettings
    import CBT.Hca.Settings.get_theory_sets;
    sets.theoryGen = get_theory_sets(sets.theoryGen); %
end

  sets.theoryGen.meanBpExt_nm = meanBpExt_nm;
   sets.theoryGen.isLinearTF = isLinearTF;
   sets.theoryGen.pixelWidth_nm = pixelWidth_nm;
  
theories = sets.theoryFold;
theorynames = sets.theoryNames;


% make theoryData folder
mkdir(sets.resultsDir);
mkdir(sets.resultsDir,timestamp);

% file where the names of theories generated in this session will be saved
matFilepathShort = fullfile(sets.resultsDir, strcat(['theories_' sprintf('%s_%s', timestamp) '.txt']));
fd = fopen(matFilepathShort,'w');    
fclose(fd);

matFastapathShort = fullfile(sets.resultsDir, strcat(['fastas_' sprintf('%s_%s', timestamp) '.txt']));
fd = fopen(matFastapathShort,'w');    
fclose(fd);

% compute free concentrations
import CBT.Hca.Core.Theory.compute_free_conc;
sets = compute_free_conc(sets);

theoryGen = struct();
tempNames = cell(1,length(sets.theoryNames));

% sets.theoryGen.meanBpExt_nm = meanBpExt_nm;
bpNmV = sets.theoryGen.meanBpExt_nm/sets.theoryGen.psfSigmaWidth_nm;

meanBpExt_nm = sets.theoryGen.meanBpExt_nm;
pixelWidth_nm = sets.theoryGen.pixelWidth_nm;
psfSigmaWidth_nm = sets.theoryGen.psfSigmaWidth_nm;
linear = sets.theoryGen.isLinearTF;
resultsDir = sets.resultsDir;

% loop over theory file folder
parfor idx = 1:length(sets.theoryNames)

%     addpath(genpath(sets.theoryFold{idx}))
    disp(strcat(['loaded theory sequence ' sets.theoryNames{idx}] ));

    % new way to generate theory, check theory_test.m to check how it works
    import CBT.Hca.Core.Theory.compute_theory_barcode;
    [theorySeq, header,bitmask] = compute_theory_barcode(sets.theoryNames{idx},sets);

	theoryBarcodes{idx} = theorySeq;
    theoryBitmasks{idx} = bitmask;

    theoryNames{idx} = header;
    theoryIdx{idx} = idx;
    bpNm{idx} = bpNmV;


    barcodeGen{idx}.rawBarcode = theorySeq;
    barcodeGen{idx}.rawBitmask = true(1,length(theorySeq));
     
        
    if sets.savetxts && ~isempty(bitmask)
        % save current theory in txt file
        C = strsplit(header(2:end),' ');
        tempNames{idx} = strcat(['theory_' C{1} '_' num2str(length(bitmask)) '_' num2str(meanBpExt_nm) '_' num2str(pixelWidth_nm) '_' num2str(psfSigmaWidth_nm) '_' num2str(linear) '_bitmask.txt']);
%         matFilename2 = strcat(['theoryTimeSeries_' C{1} '_' num2str(meanBpExt_nm) '_bpnm_barcode.txt']);
        matFilepath = fullfile(resultsDir, timestamp, tempNames{idx});
        fd = fopen(matFilepath,'w');
        fprintf(fd,' %5d', bitmask);
        fclose(fd);
    end
    
    
    
    if sets.savetxts && ~isempty(theorySeq)
        % save current theory in txt file
        C = strsplit(header(2:end),' ');
        tempNames{idx} = strcat(['theory_' C{1} '_' num2str(length(theorySeq)) '_' num2str(meanBpExt_nm) '_' num2str(pixelWidth_nm) '_' num2str(psfSigmaWidth_nm) '_' num2str(linear) '_barcode.txt']);
%         matFilename2 = strcat(['theoryTimeSeries_' C{1} '_' num2str(meanBpExt_nm) '_bpnm_barcode.txt']);
        matFilepath = fullfile(resultsDir, timestamp, tempNames{idx});
        fd = fopen(matFilepath,'w');
        fprintf(fd, strcat([' %5.' num2str(sets.theoryGen.precision) 'f ']), theorySeq);
        fclose(fd);
    end
end


for idx = 1:length(sets.theoryNames)
    if sets.savetxts && ~isempty( barcodeGen{idx}.rawBarcode)
        fd = fopen(matFilepathShort,'a'); fprintf(fd, '%s \n',fullfile(resultsDir,timestamp, tempNames{idx})); fclose(fd);
        fd = fopen(matFastapathShort,'a'); fprintf(fd, '%s \n',fullfile(theories{idx},theorynames{idx})); fclose(fd);
    end
end


% save sets
theoryGen.theoryBarcodes = theoryBarcodes;
theoryGen.theoryBitmasks = theoryBitmasks;

theoryGen.theoryNames = theoryNames;
theoryGen.theoryIdx = theoryIdx;
theoryGen.bpNm = bpNm;

theoryGen.sets = sets.theoryGen;

try
    addpath(genpath(fullfile(sets.resultsDir,[])));
catch
end
sets.theories = matFilepathShort;
sets.theory.askfortheory = 0;
sets.theory.askfornmbp = 0;
sets.theory.nmbp = sets.theoryGen.meanBpExt_nm;
sets.theory.precision = 5;
sets.theory.askfortheory = 0;


% Now compare vs theory
% get user theory
% import CBT.Hca.Settings.get_user_theory;
% [theoryStruct, sets] = get_user_theory(sets);

          % add info about barcodes as struct
    theoryStruct = cell2struct([theoryBarcodes;...
    theoryBitmasks;arrayfun(@(x) linear,1:length(theoryBarcodes),'un',false);...
    theoryNames;...
    cellfun(@(x) length(x),theoryBarcodes,'un',false);...
    arrayfun(@(x) meanBpExt_nm,1:length(theoryBarcodes),'un',false);...
    arrayfun(@(x) pixelWidth_nm,1:length(theoryBarcodes),'un',false);...
    arrayfun(@(x) psfSigmaWidth_nm,1:length(theoryBarcodes),'un',false)]',...
    {'rawBarcode','rawBitmask', 'isLinearTF','name','length','meanBpExt_nm','pixelWidth_nm','psfSigmaWidth_nm'},2);


import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct,sets );

end

