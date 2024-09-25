function theoryStructRev = gen_save_theoretical(nmPerPx, psffac, fastaFile, nmbp)
        nmpx = nmPerPx*psffac;
%         import Thry.gen_theoretical;
        [theoryStruct] = Thry.gen_theoretical(fastaFile,nmbp,0,nmpx);

        theoryStructRev = cell(1,length(theoryStruct));
        %     theoryStructF = cell(1,length(theoryStruct));

        %     datavals = [];
        %     for i=1:length(theoryStruct)
        %         datavals{i} = importdata(theoryStruct{i}.filename);
        %     end
        % add also reverse barcodes
        %     i = 1;
        vec = [];
        for i=1:length(theoryStruct)
            theoryStructRev{i}.filename =     strrep(theoryStruct{i}.filename,'barcode','rev_barcode');
            data = importdata(theoryStruct{i}.filename);
            flippedData = fliplr(data);

            %         fd = fopen( theoryStructRev{i}.filename,'w');
            %         fprintf(fd, strcat([' %5.' num2str(10) 'f ']), [data nan flippedData]);
            vec = [vec nan data nan flippedData];
            %         fclose(fd);
            %                 theoryStructF{i}.filename =     strrep(theoryStruct{i}.filename,'barcode','f_barcode');
            %
            %         fd = fopen( theoryStructF{i}.filename,'w');
            %         fprintf(fd, strcat([' %5.' num2str(10) 'f ']), [data nan(1,length(data))]);
            %         fclose(fd);

        end
        i=1;
        fd = fopen( theoryStructRev{i}.filename,'w');
        fprintf(fd, strcat([' %5.' num2str(10) 'f ']), vec);
        fclose(fd);
end