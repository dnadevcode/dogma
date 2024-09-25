function [barMat,shiftVec,barM,barMzscored] = bargroup_map_simple(bars, reFac, orBars,barcodeGen,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers)
    % Create a bargroup map from a single bargroup
    
    % Differs from previous as we use the struct to simplify the notation
    % root barcode is the first one
%     bars{1}
    
    % all the barcodes in the bargroup
    path = cellfun(@(x) str2num(x),bars);
    
    

    % check re-scaling parameters
    rescaleFactorVec = ones(1,length(path));
    rescaleFactorVec(2:end) = reFac;
    orSignVec = ones(1,length(path));
    orSignVec(2:end) = orBars;
    
    % all barcodes re-scaled with respect to first one, but we can instead
    % get mean re-scale factor which would be closer to the true factor
    meanRescaleFactor = mean(rescaleFactorVec);
    
    barS = struct();

    for ix=1:length(path)
        bar1idx=path(ix);

        barS(ix).rawBarcode = barcodeGen{bar1idx}.rawBarcode;
        barS(ix).rawBitmask = barcodeGen{bar1idx}.rawBitmask;

        if orSignVec(ix)==-1
            barS(ix).rawBarcode = fliplr(barS(ix).rawBarcode);
            barS(ix).rawBitmask = fliplr(barS(ix).rawBitmask);

        end

        rescaleFactor = rescaleFactorVec(ix)*meanRescaleFactor;
        barS(ix).rawBarcode = imresize(barS(ix).rawBarcode,'Scale' ,[1 rescaleFactor]);
        barS(ix).rawBitmask = imresize(barS(ix).rawBitmask,'Scale' ,[1 rescaleFactor]);

    %     barStr{ix} = barS;
    end
    


    % rescaleF = 0.9:0.01:1.1;
    rescaleF = 1;
    % have to save separately..
    foldSynth= strcat('output',timestamp);
    % have to save separately..
    [namesNeighbors, stI] = Core.save_bars_rescaled_txt(barS,rescaleF,foldSynth,0);
    

        %%
%     tic
    import Core.calc_overlap_first;
    [mpInbr,mp1nbr,maxMPnbr] = calc_overlap_first(barS,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,namesNeighbors);
%     toc
    maxMPnbr
    
    
    %
    strFac = stI(2:end);
    
    barIdxst = cell(1,length(stI)-1);
    for i=1:length(strFac)

        barIdxst{i} = 1*ones(1,sum(barS(1).rawBitmask));
    end
% 
%     for i=1:length(strFac)
%         A = sum(barS(i).rawBitmask);
%     %     barIdxst{i}(1:A) = (i)*ones;
%     end

    % % skip one barcode
    % [names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);

    %% calculate MP for neighbors

    %%
    for k=1:length(mp1nbr)
        [PKSnbr{k}, LOCSnbr{k}, pksUniquenbr{k}, pksUniquePosnbr{k}, bar1numnbr{k}, resfac2idxnbr{k}, or2signnbr{k}]= ...
        Core.unique_peak_barcodes(mp1nbr{k},mpInbr{k},barIdxst,strFac,k);
    end
    
    
%     import Core.mp_to_full_overlap;
%     [PCC_OVERLAPnbr, len1nbr, len2nbr,lenOverlapnbr,PCC_MPnbr] = ...
%         mp_to_full_overlap(stI(2:end),id1,id2,LOCS1nbr,mpInbr,barIdxst, MIN_OVERLAP_PIXELS, barS,1);

    import Core.create_barmat_first
    [barMat,shiftVec,barM,barMzscored] = create_barmat_first(barS,strFac,mpInbr,LOCSnbr,barIdxst);


%      figure,plot(barMzscored')
% 
%      figure,plot(barM','black')
%  
% 
% 
% bar1idx = 1;
% import Core.plot_best_match;
% [f,pos1, bar1, pos2, bar2,pcc1,pcc2,rescaleFactor, orSign] = plot_best_match(1,strFac(bar1idx),mpInbr{bar1idx},LOCSnbr{bar1idx},...
%     pksUniquePosnbr{bar1idx}, 1, barIdxst{bar1idx},1,barS,MIN_OVERLAP_PIXELS,barS(bar1idx+1));





end

