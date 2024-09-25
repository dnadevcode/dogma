function [data] = load_kymograph(tifsName,lambdasName,sets,outfold)
    %   function for loading kymographs from a structured folder
    %
    %   Args:
    %       tifsName - txt file with txt names
    %       lambdasName - txt file with lambda names
    %       sets - settings
    %       outfold - output folder
    
    if nargin < 4
        outfold = sets.outfold;
    end
    
    import CBT.Hca.Core.gen_barcodes;
    import CBT.Hca.Core.gen_consensus;
    import CBT.Hca.UI.Helper.select_all_consensus;

    % data folder
    data = cell(1,length(tifsName));
    
     for i=1:length(tifsName)
        [consensusStructs,consensusStruct, kymoStructs,barcodeGen] = create_consensus( tifsName{i}, sets.txt, sets.timeFramesNr );
        
%         [consensusStructs,consensusStruct, kymoStructs,barcodeGen] = create_consensus( tifsName{i}, sets.txt, 50 );
%         consensusStructs.treeStruct
        % consensus index, has to be more than consensus thresh
        idxConsensus = find(consensusStructs.treeStruct.maxCorCoef < sets.consensus.threshold==1, 1,'first')-1;
        numBars = cellfun(@(x) length(x.indices),consensusStruct);
        [num,idx] = max(numBars(1:idxConsensus));
        if isempty(num)
            idx = length(consensusStructs.treeStruct.maxCorCoef);
            warning("Consensus didn't pass the threshold..")
        end

% 
%         consensusStructs.consensusStruct = consensusStruct;
% 
% %         sets.output = sets.output.matDirpath;
% 
% %         if ~exists(sets.output)
% %         mkdir(sets.output);
% %         end
% 
%         import CBT.Hca.Export.export_cbc_compatible_consensus;
%         export_cbc_compatible_consensus(consensusStructs, barcodeGen,kymoStructs,sets);



        
        % run createconsensus instead? - this in case were  barcodes of
        % strange length included in the analysis, thus skewing the nm/bp
        % comparison
%         import CBT.Hca.Import.import_hca_settings;
%         [sets2] = import_hca_settings(sets.txt);
%         kymoStructsNew = kymoStructs(consensusStruct{idx}.indices);
%         barcodeGenNew =  CBT.Hca.Core.gen_barcodes(kymoStructsNew, sets2);
% %         save_kymos(barcodeGen,datetime,sets.kymosets.savekymos);
%         consensusStructsNew = CBT.Hca.Core.gen_consensus(barcodeGenNew,sets2);
% 
%         % create all consensus in hierarchical cluster
%         [consensusStructNew,~] = select_all_consensus(consensusStructsNew,sets2);

    
        % could redo consensus just with these barcodes for more accuracy..
        
        data{i}.fname = tifsName{i};
        data{i}.barcode = consensusStruct{idx}.rawBarcode;
        data{i}.bitmask = consensusStruct{idx}.rawBitmask;
        data{i}.circ = isempty(find(data{i}.bitmask==0));
        data{i}.consensusStructs = consensusStructs;
        data{i}.consensusStruct =  consensusStruct;
        data{i}.barcodeGen = barcodeGen; 
        data{i}.kymoStructs = kymoStructs; 

        % this computes theory as well, so important to check that camera
        % res is set correctly, this is sets.theoryGen.pixelWidth_nm
        try
            [lambdaValues,lengths,bppx,bpnm,bpnmTheory] = generate_lambda_values( lambdasName{i}, sets.txt,50);
             data{i}.bpnm = bpnm;
             data{i}.lengths = lengths;
            data{i}.bpnmTheory = bpnmTheory;
        catch
            data{i}.bpnm = nan;
        end
     end
    % take the first sample
    % sample1 = {'P11K0','P11K22'}; 
%     sample1 = {'P6K0','P6K25'}; 
% 
%     settings = 'sample1_settings.txt';

    % camera res - unknown (. nm/bp stretch for each sample - unknown. Hence
    % convertion px to bp - unknown. Henceforth need stretching to get the best!

%     % load data
%     data = cell(1,length(uniqueNames));
%     for i=1:length(uniqueNames)
%         data{i} = load(fullfile(uniqueFold{i},uniqueNames{i}));
%         
%            % load first
%         fname = fullfile(outfold,strcat([ uniqueNames{i} '.txt']));
%         fileID = fopen(fname,'w');
%         fprintf(fileID,'%2.16f ',data{i}.clusterConsensusData.barcode);
%         fclose(fileID);
%         data{i}.fname = fname;
%         data{i}.barcode = data{i}.clusterConsensusData.barcode;
%         data{i}.bitmask = data{i}.clusterConsensusData.bitmask;
% 
%      
%          data{i}.circ = isempty(find(data{i}.clusterConsensusData.bitmask==0));
% %         iscirculard = isempty(find(data{end}.clusterConsensusData.bitmask==0));
% 
%         % set if this comparison is circular
% %         sets.circ = iscircularq;
% 
%         % should we resize to larger size so that rescaling would have smaller affect?
%         %         data{i}.barcode = data{i}.clusterConsensusData.barcode;
%         %         data{i}.bitmask = data{i}.clusterConsensusData.bitmask;
%     end

    
end

