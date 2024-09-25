function [rescaledTable] = update_sf_barset(tableRepresentation, newSF, newO,selBarId)
    % Updates table representation with new length re-scaling factor

    %   Args:
    %       tableRepresentation - table
    %       newSF - length re-scaling factor
    %       newO - new orientation

    %   Returns:
    %       rescaledTable

    if nargin <3
        newO = 1;
    end

    if nargin < 4
        selBarId = 1;
    end

    rescaledTable = tableRepresentation;

    % only if newSF is not equal to 1
    shift = rescaledTable(selBarId,1);
    rescaledTable(1:2) = rescaledTable(1:2) - shift; % shift the table so that we re-scaling things centering on the best barcode
    
    newLens = rescaledTable(:,2)-rescaledTable(:,1)+1;

    rescaledTable(:,1) = round(rescaledTable(:,1)*newSF);
    rescaledTable(:,4) = rescaledTable(:,4)*newSF;
    rescaledTable(:,2) = round(newLens*newSF+rescaledTable(:,1)-1);

    rescaledTable(1:2) = rescaledTable(1:2) + shift;

%     rescaledTable(:,1) = rescaledTable(:,1)*newSF;
%     rescaledTable(:,4) = rescaledTable(:,4)*newSF;
%     rescaledTable(:,2) = newLens*newSF+rescaledTable(:,1)-1;

    if newO == 2
        rescaledTable(:,3) =  3-rescaledTable(:,3);

        temp1 =  -rescaledTable(:,2) ;
        temp2 = rescaledTable(:,2)-rescaledTable(:,1);
         
        rescaledTable(:,1) = temp1;
        rescaledTable(:,2) = temp1 + temp2;
    end

end

