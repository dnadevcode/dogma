function [ barcode ] = filter_barcode( barcode, sets )
    % filter_barcode
    % generates filtered kymo from kymo
    %
    %     Args:
    %         barcode: barcode to filter
    %         sets: settings, which tell as whether to filter and how
    % 
    %     Returns:
    %         barcode: filtered/original barcode
    
    % if we have chosen to filter before stretching
    % here filter size was chosen beforehand
    if sets.filter == 1
        if sets.filterMethod == 1
            barcode = barcode-imgaussfilt(barcode,sets.filterSize); % smoothening
        end
        
       if sets.filterMethod == 0
            barcode = imgaussfilt(barcode, sets.filterSize);
        end
    end
    
end

