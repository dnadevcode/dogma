function [ theoryStruct ] = convert_nmpsf_ratio( newNmBp, newPsf, theoryStruct,sets)
    % this function converts one nm/bp ratio (standard is 0.3) to another
    % (for example 0.2).
    
    % This uses the fact that convolution fo two Gaussians is a Gaussian.
    % See http://mathworld.wolfram.com/Convolution.html for details,
    % i.e. sigma_1^2+x^2 = sigma^2. So, knowing sigma_1 and sigma_2,
    % we can easily compute the unknown width x, and then just convolve
    % with it
    
    % calculate new psf
%     theoryStruct{1}.meanBpExt_nm
    pxSize = theoryStruct{1}.meanBpExt_nm/newNmBp;
    sigma1 = (1/pxSize)*theoryStruct{1}.psfSigmaWidth_nm/theoryStruct{1}.pixelWidth_nm;
    sigma =  newPsf/theoryStruct{1}.pixelWidth_nm;

    if sigma <= (sigma1-0.001)
        disp('Nm/bp ratio not converted. At least one barcode has nmbp lower than user-input' )
        return
    end
    precision = sets.theory.precision;
    matDirpath = sets.output.matDirpath;
    
    try        mkdir(fullfile(matDirpath,'theories'));
    catch
    end
        
    for i=1:length(theoryStruct)

        % first change nm to bp ratio
        fileID = fopen(theoryStruct{i}.filename,'r');
        formatSpec = '%f';
        seq = fscanf(fileID,formatSpec);
        fclose(fileID);

        % first convert to the correct length
        pxSize = theoryStruct{i}.meanBpExt_nm/newNmBp;
        
        import CBT.Core.convert_bpRes_to_pxRes;
        % convert nm ratio bitmask
        try % only if bitmask exists. first try to load
            bitname = strrep(theoryStruct{i}.filename,'barcode.txt','bitmask.txt');
            fileID = fopen(bitname,'r');
            bitmask = fscanf(fileID,' %f');
            fclose(fileID);
            bitmask = convert_bpRes_to_pxRes(bitmask, 1/pxSize);
            [~,mi,en] =fileparts(bitname);
            bitname = fullfile(fullfile(matDirpath,'theories'),strcat([mi '_converted_to' num2str(newNmBp) '_' num2str(newPsf) en ]));
            fileID = fopen(bitname,'w');
            fprintf(fileID,' %5.5f', bitmask);
            fclose(fileID);
        end
        
        % this should be correct since the convolution of two Gaussians is
        % a Gaussian, a test of this is in the script nmbpconvertion.m
        import CBT.Core.convert_bpRes_to_pxRes;
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
        sigma =  newPsf/theoryStruct{i}.pixelWidth_nm;
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
        

        theoryStruct{i}.filename = fullfile(fullfile(matDirpath,'theories'),strcat([mi '_converted_to' num2str(newNmBp) '_' num2str(newPsf) en ]));

        fileID = fopen(theoryStruct{i}.filename,'w');
        fprintf(fileID,strcat([' %2.' num2str(precision) 'f ']), seq);
        fclose(fileID);
        theoryStruct{i}.meanBpExt_nm = newNmBp;
        theoryStruct{i}.psfSigmaWidth_nm = newPsf;
        theoryStruct{i}.length = length(seq);
        

    end
end

