function [ params] = set_pval_params( )
   % based on average nm/bp ratio    
    prompt = {'Enter long theory length (px):','Enter pixel width in nm:','Length min', 'Length max','Num random barcodes','psf'};
    title = 'Input';
    dims = [1 35];
    definput = {'5610000','130','200','250','1000','300'};
    answer = inputdlg(prompt,title,dims,definput);
    try
        params.len2 = str2num(answer{1});
        params.pixelWidth_nm = str2num(answer{2});
        params.lenMin = str2num(answer{3});
        params.lenMax = str2num(answer{4});
        params.numRnd = str2num(answer{5});
        params.psfSigmaWidth_nm = str2num(answer{6});
    catch
        error('One of the parameters entered incorrectly. Enter only integers');
    end

end

