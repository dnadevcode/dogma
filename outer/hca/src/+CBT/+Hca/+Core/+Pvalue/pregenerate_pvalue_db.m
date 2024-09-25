function [] = pregenerate_pvalue_db(sets)
    % pregenerate_pvalue_db
    
    % Pregenerates p-value database file
    %     Args:
    %         file, path, params
    % 
    %     Returns:
    % 
    %     Example:
    %
    
    % if no pvalue params chosen, promt for them. Alternatively put this
    % where we generate all setting sin the beginning
    if sets.promtforparams==1
        import CBT.Hca.Import.set_pval_params;
        [ sets ] = set_pval_params( );
    end

    fullPath = strcat([sets.fold,sets.file]);

    try addpath(genpath(sets.fold));
    catch
        warning('no session file provided, creating new session file');
        %
        rS = 'presultNew.txt';
        fid = fopen(rS,'w');
        fclose(fid);
        fullPath = strcat([pwd '/' rS]);
    end
%     
%     try
%         % load p-value structure if the file already exists
%         import CBT.Hca.Import.load_pval_struct;
%         [ vals, data ] = load_pval_struct(fullPath);
%     catch
%         error('unloadable session file. Try a different session file');
%     end
          
import CBT.Hca.Core.Pvalue.precompute_pvalue_files;
precompute_pvalue_files(fullPath, sets );

end

