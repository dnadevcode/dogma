function [hcaStruct, hcaFilepath] = load_hca_struct(hcaFilepath)
    if nargin < 1
        hcaFilepath = [];
    end
    hcaSessionStruct = [];
    import CBT.Hca.Import.try_prompt_hca_filepaths;
    [~, hcaFilepaths] = try_prompt_hca_filepaths([], false);
    if isempty(hcaFilepaths)
        return;
    end
    hcaFilepath = hcaFilepaths{1};
    if isempty(hcaFilepath)
        return;
    end
    hcaStruct = load(hcaFilepath, 'hcaSessionStruct');
    if isfield(hcaStruct, 'hcaSessionStruct')
        hcaStruct = hcaStruct.hcaSessionStruct;
    end
end