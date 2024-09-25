function [ kymo ] = filter_kymos( kymo, sets )
    % filter_kymos
    % generates filtered kymo from kymo
    %
    %     Args:
    %         kymo: kymograph to filter
    %         sets: settings, which tell as whether to filter and how
    % 
    %     Returns:
    %         kymo: filtered/original kymo
    
    % if we have chosen to filter before stretching
    % here filter size was chosen beforehand
    if sets.filter == 1
        if sets.filterMethod == 1
            % todo: could vectorize filtering of rows here?
            for j=1:size(kymo,1)
                % only do the filtering for indices which are non-nan
                indd = ~isnan(kymo(j,:));
                % use imgaussfilt function to filter with a gaussian of
                % width sets.filterSize
                kymo(j,indd) = imgaussfilt(kymo(j,indd), sets.filterSize);    
            end
        end
    end
    
end

