function [T] = export_p_vals_table( theoryStruct,comparisonStructAll, barcodeGen,matDirpath )
    % export_p_vals_table

    %   Args:
    %       theoryStruct
    %       comparisonStructAll
    %       barcodeGen
    %       matDirpath

    %   Returns:
    %       Table T;
    
    fasta = cell(1,length(theoryStruct));
    for i =1:length(theoryStruct)
        if iscell(theoryStruct)  
            locs = strfind(theoryStruct{i}.filename,'/');
            fasta{i} = theoryStruct{i}.name;
        else
%             locs = strfind(theoryStruct(i).filename,'/');
            fasta{i} = theoryStruct(i).name;
        end
    end
    %fasta = cellfun(@(x) strrep(x.filename,'.txt',''), theoryStruct,'UniformOutput',false);
    if iscell(theoryStruct)  
        thrLen = cellfun(@(x) x.length, theoryStruct);
    else
        thrLen = arrayfun(@(x) theoryStruct(x).length, 1:length(theoryStruct));
    end
    
    T = table(fasta');

    for i=1:length(barcodeGen)
        maxccoef = cell2mat(cellfun(@(x) x{i}.stoufferScore(1), comparisonStructAll,'UniformOutput',0));
        lengthPx = cell2mat(cellfun(@(x) x{i}.lengthMatch,comparisonStructAll,'UniformOutput',0));
        pos = cell2mat(cellfun(@(x) x{i}.pos(1), comparisonStructAll,'UniformOutput',0));
        stretch =  cell2mat(cellfun(@(x) x{i}.bestBarStretch,comparisonStructAll,'UniformOutput',0));
        for j=1:length(pos)
           if pos(j)<= 0
               pos(j) = pos(j)+thrLen(j);
           end
        end
        maxcc = maxccoef;
        
        [d,name,ext] = fileparts(barcodeGen{i}.name);

        N = matlab.lang.makeValidName(name);
        
        T2 = table(maxcc',lengthPx', pos',stretch' ,'VariableNames',{N,strcat(['len_'  num2str(i)]),strcat(['pos_'  num2str(i)]),strcat(['stretch_'  num2str(i)])});
        T = [T T2];
    end
    disp('Saving ccvals table');
    import CBT.Hca.Export.export_cc;
    export_cc(T, matDirpath);      


end

