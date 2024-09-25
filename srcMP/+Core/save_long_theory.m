
function [names, baridx] = save_long_theory(barcodeGenGood,name)
% skips one barcode
%     if nargin < 2
%         rescaleF = 1;
%     end
%     
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(name);
%     parfor k=1:length(barskip)
%         k
k = 1;
        list =1:length(barcodeGenGood);
        
%         if ~isinf(barskip(1))
%             list(barskip(k)) = [];
%         end
        % (logical(barcodeGenGood(x).rawBitmask))
        A = arrayfun(@(x) barcodeGenGood(x).rawBarcode,list,'un',false);

        file1 = fullfile(name,strcat(num2str(k), 'long_bar.txt'));
        names{k} = fullfile(name,strcat(num2str(k), 'long_bar.txt'));

        fd = fopen(file1,'w');

        baridx{k} = nan(1,length(cell2mat(A))+length(A)-1);
        ii=1;
        for i=1:length(list)
            fprintf(fd,'%5.16f ', A{i});
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

% end

