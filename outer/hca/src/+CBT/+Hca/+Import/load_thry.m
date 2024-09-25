function [ theoryStruct,sets ] = load_thry(path, newNmBp, sets )

    if nargin < 3
                import CBT.Hca.Import.set_default_settings_code;
        sets = set_default_settings_code();
    end
    
    if nargin < 2
        path = 'allbarcodes_2018-09-05_11_31_46_session.mat';
        
        newNmBp = 0.23;
        

    end


    %%
    import CBT.Hca.UI.Helper.load_theory;
    %addpath(genpath(path))
    theoryStruct =[];
    theoryStruct = load_theory( path,theoryStruct );
   
    % compare experiment to theory
    % take the average of nm to bp from experiments, i.e. from exp file T

    import CBT.Hca.Core.Analysis.convert_nm_ratio;
    theoryStruct = convert_nm_ratio(newNmBp,theoryStruct,sets );
    sets.barcodeGenSettings.meanBpExt_nm = newNmBp;



end

