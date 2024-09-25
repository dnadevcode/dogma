function [theoryStructRev,theoryStruct,barcodeGen] = prep_thry_for_local_comp(fastaFile, nmbp, nmPerPx, psffac, islinear)
    
    % Prepares theory for local comparison
    %
    %   Args:
    %       fastaFile, nmbp, nmPerPx, psffac
    %   Returns:
    %       theoryStructRev
    
    
    if nargin < 5
        islinear = 0;
    end
    
    %
    %
        % generate theoretical
        nmpx = nmPerPx*psffac;
%         [t,matFilepathShort,theoryStruct, sets,theoryGen] = HCA_om_theory_parallel(1,0.25,sets);

        import Thry.gen_theoretical;
        [theoryStruct,~,barcodeGen] = gen_theoretical(fastaFile,nmbp,islinear,nmpx);

        theoryStructRev = cell(1,length(theoryStruct));

        timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
        [~,~] = mkdir('output');
        vec = [];
        for i=1:length(theoryStruct)
            theoryStructRev{i}.filename = fullfile('output',['rev_barcode_',timestamp,'.txt']);
            vec = [vec theoryStruct(i).rawBarcode nan fliplr(theoryStruct(i).rawBarcode) nan];
        end
        i=1;
        fd = fopen( theoryStructRev{i}.filename,'w');
        fprintf(fd, strcat([' %5.' num2str(10) 'f ']), vec);
        fclose(fd);

%     %     theoryStructF = cell(1,length(theoryStruct));
% 
%     %     datavals = [];
%     %     for i=1:length(theoryStruct)
%     %         datavals{i} = importdata(theoryStruct{i}.filename);
%     %     end
%         % add also reverse barcodes
%     %     i = 1;
%         vec = [];
%         for i=1:length(theoryStruct)
%             theoryStructRev{i}.filename =     strrep(theoryStruct{i}.filename,'barcode','rev_barcode');
%             data = importdata(theoryStruct{i}.filename);
%             flippedData = fliplr(data);
% 
%     %         fd = fopen( theoryStructRev{i}.filename,'w');
%     %         fprintf(fd, strcat([' %5.' num2str(10) 'f ']), [data nan flippedData]);
%             vec = [vec data nan flippedData nan];
%     %         fclose(fd);
%     %                 theoryStructF{i}.filename =     strrep(theoryStruct{i}.filename,'barcode','f_barcode');
%     % 
%     %         fd = fopen( theoryStructF{i}.filename,'w');
%     %         fprintf(fd, strcat([' %5.' num2str(10) 'f ']), [data nan(1,length(data))]);
%     %         fclose(fd);
% 
%         end
%         i=1;
%         fd = fopen( theoryStructRev{i}.filename,'w');
%         fprintf(fd, strcat([' %5.' num2str(10) 'f ']), vec);
%         fclose(fd);
end

% end

