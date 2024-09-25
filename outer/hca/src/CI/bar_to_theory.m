function [theoryStruct] = bar_to_theory(theorBar, theorBit, isLinearTF)
    % helping function, converts theory barcode and bitmask to theory
    % structure
        
    if ~iscell(theorBar)
        theorBar = {theorBar};
        theorBit = {theorBit};
    end
    mkdir('output');
    
    precision = 5;
    
    matDirpath =pwd;
    
    theoryStruct = cell(1,length(theorBar));
    for i=1:length(theorBar)
    
        theoryStruct{i}.filename = strcat([fullfile(matDirpath,'output','theory_') num2str(i) '_barcode.txt' ]);
        theoryStruct{i}.isLinearTF =  isLinearTF ;
        theoryStruct{i}.name = num2str(i);
        
        fileID = fopen(theoryStruct{i}.filename,'w');
        fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), theorBar{i});
        lS = length(theorBar{i});

%         if theoryStruct{i}.isLinearTF % remove the extra
%             fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), theorBar{i}(101:end-100));
%             lS = length(theorBar{i}(101:end-100));
%         else
%             fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), theorBar{i});
%             lS = length(theorBar{i});
%         end
        fclose(fileID);
        theoryStruct{i}.meanBpExt_nm = 0.3;
        theoryStruct{i}.psfSigmaWidth_nm = 300;
        theoryStruct{i}.length = lS;
        theoryStruct{i}.pixelWidth_nm = 130;
        
        bitname = strrep(theoryStruct{i}.filename,'barcode.txt','bitmask.txt');
        fileID = fopen(bitname,'w');
        fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), theorBit{i});
%         if theoryStruct{i}.isLinearTF % remove the extra
%             fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), theorBit{i}(101:end-100));
%         else
%             fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), theorBit{i});
%         end
        fclose(fileID);

    end
        
end

