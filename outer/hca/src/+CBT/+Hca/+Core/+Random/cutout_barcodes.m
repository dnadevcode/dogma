function [ barcodeGenRandom ] = cutout_barcodes( barcodeGen,sets)
    % cutout_barcodes
    % Generates random cutouts
    %     Args:
    %         barcodeGen: barcode structure
    % 
    %     Returns:
    %         barcodeGen: output barcode structure
    % 
    %     Example:
    %   

    % Generate noOfCutouts random cut-outs of length cutoutSize 
    % from the set of experimental barcodes. Bitmasks are used as input, 
    % since the actual barcodes do not matter for this purpose.
    import CBT.Hca.Core.Random.generate_random_cutouts;
    [barcodeCutOutIndices , cutoutStartPos] = generate_random_cutouts(sets.random.cutoutSize , ....
               sets.random.noOfCutouts , cellfun(@(x) x.rawBitmask, barcodeGen,'UniformOutput',false));
          
    
    % Create structure for random        
    barcodeGenRandom = cell(1,sets.random.noOfCutouts);
    % bitmask size is already computed in the settings, round it down to
    % the nearest integer
    bitmaskSizeForCutouts = round(sets.bitmasking.untrustedPx );
    


    % bitmask template for all the cutout barcodes
    bitmaskTemp=ones(1,sets.random.cutoutSize); 
    bitmaskTemp(1:bitmaskSizeForCutouts)=0;
    bitmaskTemp(end-bitmaskSizeForCutouts+1:end)=0;
    
  

    % Also explicitly extract cut-out barcodes and associated bitmasks 
    % and store in the new struct
    for counter = 1:1:sets.random.noOfCutouts

       barcodeIdxTemp = barcodeCutOutIndices(counter);
       startPosTemp = cutoutStartPos(counter);
              
        % All cut-out barcodes below can be generated from the barcode indices,  
        % start positions and cutout size, so let us add this info 
        % to a new struct (which also contains the old data).
       barcodeGenRandom{counter}.barcodeIdxTemp = barcodeIdxTemp;
       barcodeGenRandom{counter}.startPosTemp = startPosTemp;
       
       % and store the barcode
       barcodeTemp=barcodeGen{barcodeIdxTemp}.rawBarcode;
       barcodeGenRandom{counter}.rawBarcode = barcodeTemp(startPosTemp:startPosTemp+sets.random.cutoutSize-1); 
       % store the bitmask
       barcodeGenRandom{counter}.rawBitmask = bitmaskTemp;
       
       % if name already defined, add another index, which is just the idx
       % of counter
       name = strcat(['bar_' num2str(barcodeGenRandom{counter}.barcodeIdxTemp ) '_pos_' num2str(barcodeGenRandom{counter}.startPosTemp) ] );
       if ismember(name,cellfun(@(x) x.name, barcodeGenRandom(1:counter-1),'UniformOutput',false))
           name = strcat([name '_' num2str(counter)]);
       end
       barcodeGenRandom{counter}.name = name;
    end
   

%     % Plot results. Let us use the new struct to fetch the relevant info.
%     figure()
%     cutoutSize=hcaSessionCutoutStruct.cutoutSize;
%     for counter=1:1:length(barcodeGenRandom)
% 
%        barcodeIdx = barcodeGenRandom{counter}.barcodeIdxTemp;
%        plot(hcaSessionCutoutStruct.rawBarcodesCutOut{barcodeIdx}); hold on;
% 
%     end
%     hold off;
%     xlabel('pixels');
%     ylabel('intensity');

end

