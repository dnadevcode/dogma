function [theoryStructRev,theoryStruct,barcodeGen] = prep_thry_for_local_comp_simulated(lenSeq, nmbp, nmPerPx, psffac)
    
    % Prepares theory for local comparison
    %
    %   Args:
    %       lenSeq - length sequence
    % %     nmbp, nmPerPx, psffac
    %   Returns:
    %       theoryStructRev
    
    % currently only generates single sequence
    

        % generate theoretical
        nmpx = nmPerPx*psffac;
        bpPerpx = nmpx/nmbp;
        PSF_WIDTH_PIXELS = 300/nmpx;

        % rand circular barcode
        thrySeq = imgaussfilt(normrnd(0,1,1,round(lenSeq/bpPerpx)),PSF_WIDTH_PIXELS,'Padding','circular');

        theoryStructRev = cell(1,1);
        [~,~] = mkdir('barcoli');
        theoryStructRev{1}.filename = 'barcoli/bartest.txt';

        vec = [];
        for i=1:length(theoryStructRev)
            flippedData = fliplr(thrySeq);
            vec = [vec nan thrySeq nan flippedData];
        end
        % this will be used by local alignment algorithm (MP)
        fd = fopen( theoryStructRev{1}.filename,'w');
        fprintf(fd, strcat([' %5.' num2str(10) 'f ']), vec);
        fclose(fd);

        theoryStruct =[];
        barcodeGen{1}.rawBarcode = thrySeq;
        barcodeGen{1}.rawBitmask = ones(1,length(thrySeq));

end

% end

