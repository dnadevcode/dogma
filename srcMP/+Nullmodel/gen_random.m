function [fragments,barS] = gen_random(nF, psfPx, rlMin, gapPx, stdN )
    % gen_random generates random linear barcode fragments

    %   Args:
    %       nF - number fragments
    %       psfPx - point spread function in pixels
    %       rlMin - minimum (default) length
    %       gapPx - length difference between concequtive barcodes
    %       stdN - std for noise for random barcode

    %   Returns:
    %       barS - barcode cell structure with rawbarcode/rawbitmask
    %       fragments - barcode cell

    
    if nargin < 1
        nF = 100;
        psfPx = 300/50;
        rlMin = 500;
        gapPx = 10; % we generate a range of barcodes, 
    end

    if nargin < 5
        stdN = 3;
    end

    if isempty(psfPx)
        numExtra = 6;
    else
        numExtra = round(6*psfPx);
    end

    fragments = cell(1,nF);
    
    barS = cell(1,nF);

    parfor i = 1:nF
        % barlong
        barLong = stdN*randn(1, rlMin + gapPx*(i-1) + 2*numExtra);
        if ~isempty(psfPx)
            barLong = imgaussfilt(barLong,psfPx,'Padding','circular');
        end
        
        fragments{i} = barLong(numExtra+1:end-numExtra);
        
        barS{i}.rawBarcode =   fragments{i};
        barS{i}.rawBitmask =   true(size( fragments{i}));
    
    end
