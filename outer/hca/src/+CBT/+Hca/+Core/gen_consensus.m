function [ consensusStructs] = gen_consensus(barcodeGen, sets )
    % gen_consensus
    % Generates consensus for barcodes of the same length 
    %     Args:
    %         barcodeGen: barcode structure
    %         sets: Input settings to the method
    % 
    %     Returns:
    %         consensusStructs: Consensus tree structure
    % 
    %     Example:
    %          
    tic
    if sets.genConsensus == 0 || length(barcodeGen)<=1
        consensusStructs = [];
        disp('Skipping consensus generation');
        return;
    end
    
    % if doesn't exist stretchedBarcode (sets.genConsensus) // future: add
    % clustering before
    if  ~isfield(barcodeGen{1},'stretchedBarcode')
        allLengths = cellfun(@(x) length(x.rawBarcode),barcodeGen);
        commonLength = ceil(mean(allLengths));
        stretchings = commonLength./allLengths;
        strMin = min(stretchings);
        strMax =  max(stretchings);
        disp(strcat(['Barcodes are being stretched between ' num2str(strMin) ' and ' num2str(strMax)]));
        commonLength = ceil(commonLength);
        
        import CBT.Consensus.Core.convert_barcodes_to_common_length;
        import CBT.Consensus.Core.convert_bitmasks_to_common_length;

        for i=1:length(barcodeGen) % change this to a simpler function
            [barcodeGen{i}.stretchedBarcode] = cell2mat(convert_barcodes_to_common_length({barcodeGen{i}.rawBarcode}, commonLength));
            [barcodeGen{i}.stretchedrawBitmask] = cell2mat(convert_bitmasks_to_common_length({barcodeGen{i}.rawBitmask}, commonLength));
            barcodeGen{i}.stretchFactor = stretchings(i);
        end  
    end
    
        % generate consensus in case it is not aborted
	import CBT.Hca.Core.generate_consensus;
    % consensus structs: for each cluster separately 
    if isfield(barcodeGen{1},'lC')
        clusters = cellfun(@(x) x.lC,barcodeGen);
        for i=1:max(clusters)
             consensusStructs{i}  = generate_consensus( cellfun(@(x) x.stretchedBarcode,barcodeGen(clusters==i),'UniformOutput',false),...
             cellfun(@(x) x.stretchedrawBitmask,barcodeGen(clusters==i),'UniformOutput',false),...
             cellfun(@(x) x.bgMeanApprox, barcodeGen(clusters==i),'UniformOutput',false), sets );
        end
    else        

        consensusStructs  = generate_consensus( cellfun(@(x) x.stretchedBarcode,barcodeGen,'UniformOutput',false),...
             cellfun(@(x) x.stretchedrawBitmask,barcodeGen,'UniformOutput',false),...
             cellfun(@(x) x.bgMeanApprox,barcodeGen,'UniformOutput',false), sets );
    end
    
%     if length(consensusStructs)==1
%         consensusStructs = consensusStructs{1};
%     end


    timePassed = toc;
    display(strcat(['All consensuses generated in ' num2str(timePassed) ' seconds']));

end

