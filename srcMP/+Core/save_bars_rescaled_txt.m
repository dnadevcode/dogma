function [names, stridx,barcodeGenGood] = save_bars_rescaled_txt(barcodeGenGood,rescaleF,name,rev)
if nargin < 3
    name = 'bars';
end

if nargin < 4
    rev = 1;
end

if ~isstruct(barcodeGenGood)
barcodeGenGood = cell2struct([cellfun(@(x) double(x.rawBarcode),barcodeGenGood,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGenGood,'un',false)]',{'rawBarcode','rawBitmask'},2);
end

[SUCCESS,MESSAGE,MESSAGEID] = mkdir(name);

stridx = cell(1,length(barcodeGenGood));
names = cell(1,length(barcodeGenGood));

for i=1:length(barcodeGenGood)
    names{i} = fullfile(name,strcat(num2str(i), '_bar.txt'));

    file1= fullfile(name,strcat(num2str(i), '_bar.txt'));

    A = arrayfun(@(y) imresize(barcodeGenGood(i).rawBarcode(logical(barcodeGenGood(i).rawBitmask)),'Scale' ,[1 y]),rescaleF,'un',false);
    
    fd = fopen(file1,'w');
    
    ii = 1;
    if rev==1
        stridx{i} = nan(1,2*(sum(cellfun(@(y) length(y),A)) + length(A)-1));
    else
        stridx{i} = nan(1,(sum(cellfun(@(y) length(y),A)) + length(A)-1));
    end

    for j=1:length(A)
        fprintf(fd,'%5.16f ', A{j});
        fprintf(fd,'%5.16f ', NaN);
        stridx{i}(ii:ii+length(A{j})-1) = j;
        ii = ii+length(A{j})+1;
    end
    
    if rev==1
        for j=1:length(A)
            fprintf(fd,'%5.16f ', fliplr(A{j}));
            fprintf(fd,'%5.16f ', NaN);
            stridx{i}(ii:ii+length(A{j})-1) = -j;
            ii = ii+length(A{j})+1;
        end
    end
    fclose(fd);
   
end


end

