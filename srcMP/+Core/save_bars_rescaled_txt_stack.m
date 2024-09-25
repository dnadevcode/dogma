function [names, stridx, baridx] = save_bars_rescaled_txt_stack(barcodeGenGood,rescaleF,name,rev)
if nargin < 3
    name = 'bars';
end

if nargin < 4
    rev = 1;
end
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(name);

names{1} = fullfile(name,'rescaled_bar.txt');

file1= fullfile(name,'rescaled_bar.txt');
    
stridx = cell(1,length(barcodeGenGood));

fd = fopen(file1,'w');

for i=1:length(barcodeGenGood)
    A = arrayfun(@(y) imresize(barcodeGenGood(i).rawBarcode(logical(barcodeGenGood(i).rawBitmask)),'Scale' ,[1 y]),rescaleF,'un',false);
    
    ii = 1;
    if rev==1
        stridx{i} = nan(1,2*(sum(cellfun(@(y) length(y),A)) + length(A)));
        baridx{i} = i*ones(1,2*(sum(cellfun(@(y) length(y),A)) + length(A)));
    else
        stridx{i} = nan(1,(sum(cellfun(@(y) length(y),A)) + length(A)));
        baridx{i} = i*ones(1,(sum(cellfun(@(y) length(y),A)) + length(A)));
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
end
    fclose(fd);

    stridx = cell2mat(stridx);
    baridx = cell2mat(baridx);
    


end

