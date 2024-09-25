function [theoryStruct] = load_theory_into_txts(theoryStruct,theoryFold, bars,bits, names, meanbpnm, pixelWidth_nm, psfSigmaWidth_nm, isLinearTF, precision)

    for i=1:length(bars)
        fname = fullfile(theoryFold, strcat(num2str(i),'_barcode.txt'));
        fileID = fopen(fname,'w');
        % choose the precision from settings
        fprintf(fileID,strcat(['%2.' num2str(precision) 'f ']),bars{i});
        fclose(fileID);
        try
            fnameb = strrep(fname,'_barcode','_bitmask');
            fileID = fopen(fnameb,'w');
            % choose the precision from settings
            fprintf(fileID,strcat(['%2.' num2str(precision) 'f ']),bits{i});
            fclose(fileID);
        catch
        end
        theoryStructN{i}.filename = fname;
        theoryStructN{i}.name = names{i};
        theoryStructN{i}.length =  length(bars{i});
        theoryStructN{i}.meanBpExt_nm = meanbpnm;
        theoryStructN{i}.pixelWidth_nm = pixelWidth_nm;
        theoryStructN{i}.psfSigmaWidth_nm = psfSigmaWidth_nm;
        theoryStructN{i}.isLinearTF = isLinearTF;
    end
    theoryStruct = [theoryStruct theoryStructN];
end

