function [ filterSize] = get_hca_filter_settings(default_filterSize,titleText)
    
    filterSize = default_filterSize;
    options.Resize = 'on';
    % sets number of time-frames
    answer = inputdlg(...
        { ...
            'Choose filter size', .....
        }, ... % prompt
        titleText, ... % dialog title
        1, ... % number of lines
        { ...
            num2str(filterSize), ...
        },options);
    if ~isempty(answer)
        filterSize = str2double(answer{1});
    end

    isAcceptable = struct;
    isAcceptable.filterSize = (filterSize >= 0);
  
    if not(isAcceptable.filterSize)
        filterSize = default_filterSize;
        warning('Bad input for filter size! Try again!');
    end
end