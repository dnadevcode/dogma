function [overlapStruct2,mpI1,mp1,stridx,baridx2] = calc_overlap_mp_multidim(barcodeGen, bT,sF, MIN_OVERLAP_PIXELS,timestamp,scampDir,numWorkers)
    % main fun - calculated MP overlaps
%     if nargin < 5
%         if ispc % for PC get from the initialization
%             sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP' ;
%         else
%             sets.SCAMP_LINE = '~/SCAMP/';
%         end
        numWorkers = 32;
%     else
%         sets.SCAMP_LINE = uigetdir("Please provide SCAMP folder");
%     end

    SCAMP_LINE = 'SCAMP';

    % TODO: alternative, matlab's version of mp algorithm in case scamp is
    % not available.

    
    %% Same with MP
    barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
        cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
    % have to save separately..
    foldSynth=strcat('output',timestamp);
    [~,~] =mkdir(foldSynth);

    sF = 1;
    % have to save separately.. %todo: save with mask
    [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);
    % skip one barcode
    [names2, baridx2] = Core.save_long_bar_txt_skip(bT,inf,foldSynth);


        % tic
    import Core.calc_overlap;
    [mpI1, mp1, maxMP] = calc_overlap([], timestamp, SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers, namesBar, repmat(names2,1,length(barStruct)));
    % toc

    % find max along mp1:
    mps = cell2mat(mp1);
    figure,imagesc(mps')
    xlabel('Pos along theory')
    ylabel('Barcode')
    title('Local alignment score')

    mpI = cell2mat(mpI1);

    % todo:
    % each mp/mpI location convert 

    for ii=1:length(baridx2)
        baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
    end

    % we want to create a nicer structure for all-to-all comparison and
    % contains easily accesible data.
    % tic
%     import Core.mp_res_to_struct;
%     [overlapStruct2] = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
    % toc

    import Core.mp_to_struct;
    [overlapStruct2] = mp_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
  

% import Core.mp_res_to_struct_mini;
% [overlapStruct2] = mp_res_to_struct_mini(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);



end

