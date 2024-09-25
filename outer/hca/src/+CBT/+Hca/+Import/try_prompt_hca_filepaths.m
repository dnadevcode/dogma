function [aborted, hcaFilepaths] = try_prompt_hca_filepaths(promptTitle, multiselectOn)
    % This is quite a general function and should not be made every time a
    % new just by changing prompt title. So next time this is needed do not
    % create a new one.
    
    if (nargin < 1)
        promptTitle = [];
    end
    if (nargin < 2) || isempty(promptTitle)
        multiselectOn = true;
    end
    if isempty(promptTitle)
        if multiselectOn
            promptTitle = 'Select HCA session Files';
        else
            promptTitle = 'Select HCA session File';
        end
    end
    if multiselectOn
        multiselectStr = 'on';
    else
        multiselectStr = 'off';
    end
    hcaFilepaths = cell(0, 1);

    [hcaFilenames, dirpath] = uigetfile({'*.mat;'; '*.txt;'}, promptTitle, 'MultiSelect', multiselectStr);

    aborted = isequal(dirpath, 0);
    if aborted
        return;
    end
    if not(iscell(hcaFilenames))
        hcaFilenames = {hcaFilenames};
    end
    hcaFilepaths = fullfile(dirpath, hcaFilenames(:));
end