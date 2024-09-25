function [ meanBpExt_nm] = get_nmbp_settings(default_meanBpExt_nm,titleText)
    
    meanBpExt_nm = default_meanBpExt_nm;
    options.Resize = 'on';
    % sets number of time-frames
    answer = inputdlg(...
        { ...
            'Choose mn/bp extension', .....
        }, ... % prompt
        titleText, ... % dialog title
        1, ... % number of lines
        { ...
            num2str(meanBpExt_nm), ...
        },options);
    if ~isempty(answer)
        meanBpExt_nm = str2double(answer{1});
    end

    isAcceptable = struct;
    isAcceptable.meanBpExt_nm = (meanBpExt_nm >= 0);
  
    if not(isAcceptable.meanBpExt_nm)
        meanBpExt_nm = default_meanBpExt_nm;
        warning('Bad input for nm/bp extension! Try again!');
    end
end