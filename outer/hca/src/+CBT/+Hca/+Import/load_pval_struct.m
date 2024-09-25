function [ vals, data ] = load_pval_struct( rS )
    % load_pval_struct
    %
    % Loads p-value structure
    %
    %     Args:
    %         rS
    % 
    %     Returns:
    %         vals, data 
    %     Example:
    
    
    % check if file exists
    try
        fid = fopen(rS,'a');
        fclose(fid);
    catch
        fid = fopen(rS,'w');
        fclose(fid);
    end


    fid = fopen(rS);
    lines = textscan(fid, '%s', 'Delimiter','\n');
    newRes = cell(1,length(lines{1}));
    for i=1:length(lines{1})
        newRes{i} = textscan(lines{1}{i},'%f', 'Delimiter',' ');
    end
    fclose(fid);
    
    vals = cellfun(@(x) x{1}(1), newRes);
    data = cellfun(@(x) x{1}(2:end), newRes, 'UniformOutput',false);

end

