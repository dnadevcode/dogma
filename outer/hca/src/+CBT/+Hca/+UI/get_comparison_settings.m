function [ meanBpExt_nm,concNetropsin_molar,concYOYO1_molar,concDNA,psfSigmaWidth_nm,pixelWidth_nm,deltaCut,isLinearTF,widthSigmasFromMean,computeFreeConcentrations,model] = get_comparison_settings(default_meanBpExt_nm)

    import CBT.get_default_barcode_gen_settings;
    barcodeGenSettings = get_default_barcode_gen_settings();
    
    meanBpExt_nm = barcodeGenSettings.meanBpExt_nm;
    concNetropsin_molar = barcodeGenSettings.concNetropsin_molar;
    concYOYO1_molar = barcodeGenSettings.concYOYO1_molar;
    concDNA = barcodeGenSettings.concDNA;
    psfSigmaWidth_nm = barcodeGenSettings.psfSigmaWidth_nm;
    pixelWidth_nm = barcodeGenSettings.pixelWidth_nm;
    deltaCut = barcodeGenSettings.deltaCut;
    isLinearTF = barcodeGenSettings.isLinearTF;
    widthSigmasFromMean = barcodeGenSettings.widthSigmasFromMean;
    computeFreeConcentrations =  barcodeGenSettings.computeFreeConcentrations;
    model = barcodeGenSettings.model;
    % Set a threshold for grouping barcodes.
    options.Resize = 'on';
    answer = inputdlg(...
        { ...
            'Choose nm/bps','concNetropsin Mmolar','concYOYO1 Mmolar','conc DNA Mmolar', 'psfSigmaWidth nm', 'pixelWidth_nm','deltaCut','isLinearTF','widthSigmasFromMean','compute free concentrations','model' ...
        }, ... % prompt
        'Parameter choice for barcode generation', ... % dialog title
        1, ... % number of lines
        { ...
            num2str(meanBpExt_nm), num2str(concNetropsin_molar),num2str(concYOYO1_molar),num2str(concDNA),num2str(psfSigmaWidth_nm),num2str(pixelWidth_nm),num2str(deltaCut),num2str(isLinearTF),num2str(widthSigmasFromMean),num2str(computeFreeConcentrations),num2str(model)...
        }, options... % default value
    );

    if ~isempty(answer)
        meanBpExt_nm = str2double(answer{1});
        concNetropsin_molar = str2double(answer{2});
        concYOYO1_molar = str2double(answer{3});
        concDNA = str2double(answer{4});
        psfSigmaWidth_nm =  str2double(answer{5});
        pixelWidth_nm =  str2double(answer{6});
        deltaCut =  str2double(answer{7});
        isLinearTF =  str2double(answer{8});
        widthSigmasFromMean =  str2double(answer{9});
        computeFreeConcentrations = str2double(answer{10});
        model = answer{11};
    end
    
    isAcceptable = struct;
    isAcceptable.meanBpExt_nm = (meanBpExt_nm >= 0);
  

    if not(isAcceptable.meanBpExt_nm)
        meanBpExt_nm = default_meanBpExt_nm;
        warning('Bad input for nm/bp! Try again!');
    end
  
end