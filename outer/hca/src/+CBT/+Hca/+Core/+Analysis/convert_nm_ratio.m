function [ theoryStruct ] = convert_nm_ratio( newNmBp, theoryStruct,sets)
    % this function converts one nm/bp ratio (standard is 0.3) to another
    % (for example 0.2).
    
    % This uses the fact that convolution fo two Gaussians is a Gaussian.
    % See http://mathworld.wolfram.com/Convolution.html for details,
    % i.e. sigma_1^2+x^2 = sigma^2. So, knowing sigma_1 and sigma_2,
    % we can easily compute the unknown width x, and then just convolve
    % with it

    %   Args:  
    %       newNmBp - new nm/bp ratio
    %       theoryStruct - structure containing txts with theories
    %       sets - settings
    %
    %   Returns:
    %       theoryStruct - theory structure with updated names
    
    if iscell(theoryStruct)
        minNmBp = min(cellfun(@(x) x.meanBpExt_nm,theoryStruct));
    else
        minNmBp = min(arrayfun(@(x) theoryStruct(x).meanBpExt_nm,1:length(theoryStruct)));
    end
       
        
    if newNmBp >= (minNmBp-0.001)
        disp('Nm/bp ratio not converted. At least one barcode has nmbp lower than user-input' )
        return
    end


    import CBT.Core.convert_bpRes_to_pxRes;

    if sets.theory.theoryDontSaveTxts % new in 4.7.0
        % all theories SHOULD have the same meanBpExt_nm, so move these
        % outside of loop
        pxSize = theoryStruct(1).meanBpExt_nm/newNmBp;
        sigma1 = (1/pxSize)*theoryStruct(1).psfSigmaWidth_nm/theoryStruct(1).pixelWidth_nm;
        sigma =  theoryStruct(1).psfSigmaWidth_nm/theoryStruct(1).pixelWidth_nm;
        % size of final sigma
        sigmaDif = sqrt(sigma^2-sigma1^2);
       
        %
        if length(theoryStruct) > 30
            parfor i=1:length(theoryStruct)
                if theoryStruct(i).length ~= 0
                    % first convert to the correct length
                    seq = convert_bpRes_to_pxRes(theoryStruct(i).rawBarcode, 1/pxSize);
            
                    % length of kernel
                    hsize = size(seq,2);
                    
                    % kernel
                    ker = circshift(images.internal.createGaussianKernel(sigmaDif, hsize),round(hsize/2));   
                    
                    % conjugate of kernel in phase space
                    multF=conj(fft(ker'));
                    % convolved with sequence ->
                    theoryStruct(i).rawBarcode = ifft(fft(seq).*multF); 
                    
                    if isfield(theoryStruct(i),'rawBitmask')&&~isempty(theoryStruct(i).rawBitmask)
                        theoryStruct(i).rawBitmask = convert_bpRes_to_pxRes(theoryStruct(i).rawBitmask, 1/pxSize);
                    end
                    theoryStruct(i).meanBpExt_nm = newNmBp;
                    theoryStruct(i).length = hsize;
                 
    
                end
            end
        else % do not run parfor for small datasets
            for i=1:length(theoryStruct)
                if theoryStruct(i).length ~= 0
                    % first convert to the correct length
                    seq = convert_bpRes_to_pxRes(theoryStruct(i).rawBarcode, 1/pxSize);

                    % length of kernel
                    hsize = size(seq,2);

                    % kernel
                    ker = circshift(images.internal.createGaussianKernel(sigmaDif, hsize),round(hsize/2));

                    % conjugate of kernel in phase space
                    multF=conj(fft(ker'));
                    % convolved with sequence ->
                    theoryStruct(i).rawBarcode = ifft(fft(seq).*multF);

                    if isfield(theoryStruct(i),'rawBitmask')&&~isempty(theoryStruct(i).rawBitmask)
                        theoryStruct(i).rawBitmask = convert_bpRes_to_pxRes(theoryStruct(i).rawBitmask, 1/pxSize);
                    end
                    theoryStruct(i).meanBpExt_nm = newNmBp;
                    theoryStruct(i).length = hsize;


                end
            end
        end
    else % old, if need to save to individual txts (for running with tools outside hca)
        precision = sets.theory.precision;
        matDirpath = sets.output.matDirpath;
        [~,~] =mkdir(fullfile(matDirpath,'theories'));
        parfor i=1:length(theoryStruct)
            if theoryStruct{i}.length ~= 0
        
                % first change nm to bp ratio
                fileID = fopen(theoryStruct{i}.filename,'r');
                formatSpec = '%f';
                seq = fscanf(fileID,formatSpec);
                fclose(fileID);
        
                % first convert to the correct length
                pxSize = theoryStruct{i}.meanBpExt_nm/newNmBp;
                
                % convert nm ratio bitmask
                try % only if bitmask exists. first try to load
                    bitname = strrep(theoryStruct{i}.filename,'barcode.txt','bitmask.txt');
                    fileID = fopen(bitname,'r');
                    bitmask = fscanf(fileID,' %f');
                    fclose(fileID);
                    bitmask = convert_bpRes_to_pxRes(bitmask, 1/pxSize);
                    [~,mi,en] =fileparts(bitname);
                    bitname = fullfile(fullfile(matDirpath,'theories'),strcat([mi '_converted_to' num2str(newNmBp) en ]));
                    fileID = fopen(bitname,'w');
                    fprintf(fileID,' %5.5f', bitmask);
                    fclose(fileID);
                catch
                    % do nothing
                end
                
                % this should be correct since the convolution of two Gaussians is
                % a Gaussian, a test of this is in the script nmbpconvertion.m
                seq = convert_bpRes_to_pxRes(seq, 1/pxSize);
        %         seq = interp1(1:length(seq),seq,linspace(1,length(seq),(length(seq)*1/pxSize)));
                
        %         figure,plot(seq1)
        %         hold on
        %         plot(seq2)
        %         legend({'converted to bp','Interpolated'})
        %         
        %         % We assume bitmasks to be only ones for theories
                %hcaSessionStruct.theoryGen.bitmask{i} = convert_bpRes_to_pxRes(hcaSessionStruct.theoryGen.bitmask{i}, 1/pxSize);
                sigma1 = (1/pxSize)*theoryStruct{i}.psfSigmaWidth_nm/theoryStruct{i}.pixelWidth_nm;
                sigma =  theoryStruct{i}.psfSigmaWidth_nm/theoryStruct{i}.pixelWidth_nm;
                % size of final sigma
                sigmaDif = sqrt(sigma^2-sigma1^2);
               
                % length of kernel
                hsize = size(seq,2);
           
                % kernel
                ker = circshift(images.internal.createGaussianKernel(sigmaDif, hsize),round(hsize/2));   
                
                % conjugate of kernel in phase space
                multF=conj(fft(ker'));
        
                % convolved with sequence ->
                seq = ifft(fft(seq).*multF); 
                %fname = strcat(['fold/theory_' barcodeData.hcaSessionStruct.theoryNames{i} '.txt']);
                [~,mi,en] =fileparts(theoryStruct{i}.filename);
                
    %             try        
    %             catch
    %             end
                theoryStruct{i}.filename = fullfile(fullfile(matDirpath,'theories'),strcat([mi '_converted_to' num2str(newNmBp) en ]));
        
                fileID = fopen(theoryStruct{i}.filename,'w');
                fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), seq);
                fclose(fileID);
                theoryStruct{i}.meanBpExt_nm = newNmBp;
                theoryStruct{i}.length = length(seq);
                
            end
        end
    end
end

