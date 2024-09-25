function [T] = export_infoscore_vals_table( barcodeGenC,matDirpath )
    % export_infoscore_vals_table

    % Input:
       % barcodeGenC,matDirpath
    % Output:

    % barcode names
    T = table(cellfun(@(x) matlab.lang.makeValidName(x.name) ,barcodeGenC,'UniformOutput',false)','VariableNames',{'experiment'});

	T2 = table(cellfun(@(x) x.infoscore.mean,barcodeGenC)',cellfun(@(x) x.infoscore.std,barcodeGenC)' ,'VariableNames',{'mean','std'});

    
    % score1: mean+3std
    TS1 =  table(cellfun(@(x) x.infoscore.mean,barcodeGenC)'+3*cellfun(@(x) x.infoscore.std,barcodeGenC)' ,'VariableNames',{'score1'});

    
    % score2: mean/std
    TS2 =  table(cellfun(@(x) x.infoscore.mean,barcodeGenC)'./cellfun(@(x) x.infoscore.std,barcodeGenC)' ,'VariableNames',{'score2'});

    is3 = zeros(1,length(barcodeGenC));
    for i=1:length(barcodeGenC);
        vals =barcodeGenC{i}.rawBarcode(logical(barcodeGenC{i}.rawBitmask));
        is3(i) = sum(abs(zscore(vals)) > 3)/length(vals);
    end
    % score3: proportion of values outside mean+3std
    TS3 =  table(is3' ,'VariableNames',{'score3'});

    % score4:  FS = variance in signal region/variance for background.  
    %Since a barcode is really only meaningful if the
    %amplitudes in the signal regions is (much) larger than the noise
    %amplitude, we could calculate
    is4 = zeros(1,length(barcodeGenC));
    for i=1:length(barcodeGenC)
        is4(i) =std(barcodeGenC{i}.rawBarcode(logical(barcodeGenC{i}.rawBitmask)))/barcodeGenC{i}.bgStdApprox;
    end
	TS4 =  table(is4' ,'VariableNames',{'FS'});

	is5 = zeros(1,length(barcodeGenC));
    for i=1:length(barcodeGenC)
        der =diff(barcodeGenC{i}.rawBarcode(logical(barcodeGenC{i}.rawBitmask)));
        stder = std(der);
        try 
            [PKS,~] = findpeaks(der, 'MinPeakHeight',3*stder);
            is5(i) = length(PKS);
        catch
            is5(i) = nan;
        end
    end
    
	TS5 =  table(is5' ,'VariableNames',{'NOALJ'});


    T = [T T2 TS1 TS2 TS3 TS4 TS5];

    disp('Saving infoscore table');
    disp('Score1: mean+3std');
    disp('Score2: mean/std');
    disp('Score3: proportion of values outside mean+3std');
	disp('Score4: FS');
    disp('Score5: NOALJ.');


    import CBT.Hca.Export.export_cc;
    export_cc(T, matDirpath);      


end

