function [maxPCC, allPCC] = gen_scores_pcc(...
    NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,RAND_LENGTH_2,...
    NN,sF,out, gap)

    if nargin < 8
        gap = 10;
    end
        
    [~,~]=mkdir(out);
    
    maxPCC=cell(1,NN);
    
    
    comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),0,50);

    
    import Nullmodel.gen_random;

    if nargout >= 2
        allPCC = cell(1,NN);
    end
    
    pccs = cell(1,NN);
    for k=1:NN
        maxPCC{k} = nan(1,NUM_RAND_FRAGMENTS);
        
        [barcodes] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,0);
    
        %
        bitmasks = cellfun(@(x) true(size(x)), [barcodes], 'un', 0);
    
    
        % create bar structure 
        barSynth = cell2struct([barcodes; bitmasks]',{'rawBarcode','rawBitmask'},2);
    
        [barcodes2] = gen_random(1,PSF_WIDTH_PIXELS,RAND_LENGTH_2+(k-1)*gap,0);
        bitmasks2 = cellfun(@(x) true(size(x)), [barcodes2], 'un', 0);
    
        barSynth2 = cell2struct([barcodes2; bitmasks2]',{'rawBarcode','rawBitmask'},2);
    
    
        % barSynth = [cell2struct(barcodes,{'rawBarcode'},1);cell2struct(bitmasks,{'bitmasks'},1)];
        pccs{k} = zeros(1,length(barSynth));
        for i=1:length(barSynth)
            A = arrayfun(@(y) imresize(barSynth(i).rawBarcode(barSynth(i).rawBitmask),'Scale' ,[1 y]),sF,'un',false);
            pcCur = 0;
            for j=1:length(A)
                if nargout >= 2
                    [pccThis,~,~,~,~,allPCC{k}{i}{j}] = comparisonFun(A{j},barcodes2{1},ones(1,length(A{j})),ones(1,length(barcodes2{1})),2^16);
                else
                    pccThis = comparisonFun(A{j},barcodes2{1},ones(1,length(A{j})),ones(1,length(barcodes2{1})),2^16);
                end
                pcCur = max([pcCur,pccThis]);
            end
            pccs{k}(i) = pcCur;
        end
        % have to save separately..
    
        maxPCC{k}= pccs{k};
    

    end


end

