function [valRes,val,sorti] = scaled_specific_barcode_vs_theory(compI,selBarId,barcodeIslandsData,theoryStruct,idx,forceDirection)

% if nargin < 6
%     forceDirection = 0;
% end
    % forceDirection - make sure that the reference is oriented the same
    % way as the theory
import Core.update_sf_barset;
import Core.synth_to_table;




    %% extract tables
    selectedBarStretch = compI{selBarId}.bestBarStretch;  % best bar stretch (bbS) along the theory
    selectedBarOr = (compI{selBarId}.or(1)~=barcodeIslandsData{idx}(selBarId,3))+1; % orientation w.r.t. theory

    % Update EXPERIMENT TABLE barcodeIslandsData{idx} to correct scaling factor & orientation
    [expTable] = update_sf_barset(barcodeIslandsData{idx}, selectedBarStretch/barcodeIslandsData{idx}(selBarId,4), selectedBarOr, selBarId);

    % Convert all detected locations to tableS
    [tableReference] = synth_to_table(compI);

    % update reference table to have the same orientation as the
    [tableReferenceUpd] = update_sf_barset(tableReference, tableReference(selBarId,4)/selectedBarStretch, 1); % should this be for specific bar since we change scaling?

    % experimental island
%     if forceDirection && tableReferenceUpd(selBarId,3) == 2
%         [expTable] = update_sf_barset(expTable, 1, 2);
%         [tableReferenceUpd] = update_sf_barset(tableReferenceUpd, 1, 2);
%     end

    % align experimental table to match with reference table
    expTable(:,1:2) = expTable(:,1:2)-expTable(selBarId,1)+tableReferenceUpd(selBarId,1); %
   
    % positions along exp table
    bestPosFound = expTable(:,1);
    
    % positions along reference table
    posTrue = tableReferenceUpd(:,1);

    % pDif for each position (sqrt of square of differences)
    pDif{idx} =  min([sqrt((posTrue-bestPosFound).^2) sqrt((posTrue+theoryStruct.length-bestPosFound).^2)...
        sqrt((posTrue-theoryStruct.length-bestPosFound).^2)]');

    [val,sorti] = sort(posTrue);

    valRes.pDif = pDif{idx};
%     valRes{idx}.compI = compI;
    valRes.bbSAll = expTable;
%     valRes{idx}.bars = bars;
    valRes.bestPosFound = bestPosFound;
    valRes.selBarId = selBarId;
    valRes.tableSUpd = tableReferenceUpd;
    valRes.tableS = tableReference;
    valRes.posTrue = posTrue;

end

