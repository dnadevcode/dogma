function [names, baridx] = save_long_bar_txt_skip(barcodeGenGood, barskip,name)
% skips one barcode
%     if nargin < 2
%         rescaleF = 1;
%     end

    if ~isstruct(barcodeGenGood)
    barcodeGenGood = cell2struct([cellfun(@(x) double(x.rawBarcode),barcodeGenGood,'un',false);...
        cellfun(@(x) x.rawBitmask,barcodeGenGood,'un',false)]',{'rawBarcode','rawBitmask'},2);
    end

%     
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(name);
    names = cell(1,length(barskip));
    baridx = cell(1,length(barskip));
    for k=1:length(barskip)
%         k
        list =1:length(barcodeGenGood);
        
        if ~isinf(barskip(1))
            list(barskip(k)) = [];
        end

        if ~isempty(barcodeGenGood(list(1)).rawBitmask) % if at least one does not contain a bitmask, we skip that
            A = arrayfun(@(x) barcodeGenGood(x).rawBarcode(logical(barcodeGenGood(x).rawBitmask)),list,'un',false);
        else
            A = arrayfun(@(x) barcodeGenGood(x).rawBarcode,list,'un',false);
        end

        file1 = fullfile(name,strcat(num2str(k), 'long_bar.txt'));
        names{k} = fullfile(name,strcat(num2str(k), 'long_bar.txt'));

        fd = fopen(file1,'w');

        baridx{k} = nan(1,length(cell2mat(A))+length(A)-1);
        ii=1;
        for i=1:length(list)
            fprintf(fd,'%5.16f ', A{i}); % save in less precision?
            fprintf(fd,'%5.16f ', NaN);
            baridx{k}(ii:ii+length(A{i})-1) = list(i);
            ii = ii+length(A{i})+1;
        end
        fclose(fd);
%         file2= fullfile('barsLong',strcat(num2str(k), '_barF.txt'));
%         namesF{k} = fullfile('barsLong',strcat(num2str(k), '_barF.txt'));
% 
%         fd = fopen(file2,'w');
% 
%         for i=1:length(list)
%             fprintf(fd,'%5.16f ',fliplr(A{i}));
%             fprintf(fd,'%5.16f ', NaN);
%         end
%         fclose(fd);
    end

end

