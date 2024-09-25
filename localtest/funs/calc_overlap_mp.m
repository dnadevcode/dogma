function [overlapStruct2,mpI1,mp1,stridx,baridx2] = calc_overlap_mp(barcodeGen,sF, MIN_OVERLAP_PIXELS,timestamp,numWorkers)
    % Calculating overlaps between all pairs of barcodes using SCAMP matrix
    % profile algorithm
    
    %   Args:
    %       barcodeGen - barcode structure with rawBarcode and rawBitmask
    %       fields
    %       sF - length re-scaling factors
    %       MIN_OVERLAP_PIXELS - minimum overlap of pixels between any two
    %       barcodes
    %       timestamp - timestamp for temporary output folder
    %       numWorkers - number of workers to use for calculation
    
    %   Returns:
    %       overlapStruct2 - overlap structure, with fields of 
    %         score
    %         fullscore
    %         overlaplen
    %         pA
    %         pB
    %         or
    %         bestBarStretch
    %         lenA
    %         lenB
    %         partialScore
    %         partialLength
    %         h


    % main fun - calculated MP overlaps
    if nargin < 5
        numWorkers = 32;
    end

    SCAMP_LINE = 'SCAMP';

    % TODO: alternative, matlab's version of mp algorithm in case scamp is
    % not available.

    
    %% Same with MP
    if ~isstruct(barcodeGen)
        barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
            cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
    else
        barStruct = barcodeGen;
    end
   
    % have to save separately..
    foldSynth = tempname('output'); %unique folder name. to be removed after calculation
    [~,~] = mkdir(foldSynth);

    % have to save separately..
    [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);
    % skip one barcode
    [names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);
    
    % tic
    import Core.calc_overlap;
    [mpI1, mp1, maxMP] = calc_overlap(barStruct, foldSynth, SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers, names2, namesBar);
    % toc


    for ii=1:length(baridx2)
        baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
    end

    % use parallel here to speed up
    import Core.mp_to_struct;
    [overlapStruct2] = mp_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);

    [~,~] =  rmdir(foldSynth,'s'); % remove temporary folder


end

