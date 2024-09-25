function [ rezMaxM,bestBarStretch,bestLength ] = on_compare(barcodeGen,theoryStruct,comparisonMethod,stretchFactors,w,numPixelsAroundBestTheoryMask, positionsOnThry, sets)
    % on_compare_theory_to_exp
    % Compares experiments to single theory
    %     Args:
    %         sets: settings structure
    %         barcodeGen: barcode structure
    %         theoryStruct: theory structure
    % 
    %     Returns:
    %         comparisonStructure: comparison structure

    % create the rezult structure
%     comparisonStructure = cell(1, length(barcodeGen));
    %ccM = cell(length(rawBarcodes),1);

    switch comparisonMethod
        case 'spearman'
            import SignalRegistration.unmasked_spearman_corr;
            comparisonFun = @(x,y,w,z,u) unmasked_spearman_corr(x',y',w);
        case 'spearman2'
            import SignalRegistration.unmasked_spearman_corr2;
            comparisonFun = @(x,y,w,z,u) unmasked_spearman_corr2(x',y',w);
        case 'spearman3'
            import SignalRegistration.unmasked_spearman_corr3;
            comparisonFun = @(x,y,w,z,u) unmasked_spearman_corr3(x',y',w);

        case 'unmasked_pcc_corr'
            import SignalRegistration.unmasked_pcc_corr;
            comparisonFun = @(x,y,z,w,u) unmasked_pcc_corr(x,y,z);
%             import CBT.Hca.UI.Helper.get_best_parameters;

%             parameterfun = @(x)
        case 'mass_dot'
            comparisonFun = @(x,y,z,w,u) unmasked_MASS_DOT_CC(y,x,z,w,2^(4+nextpow2(length(x))),theoryStruct.isLinearTF,numPixelsAroundBestTheoryMask);            
       case 'masked_pcc'
            % choose k just higher than the length of small sequence for
            % best precision. (larger k though could increase speed)
            comparisonFun = @(x,y,z,w,u) masked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),theoryStruct.isLinearTF,numPixelsAroundBestTheoryMask);

        case 'mass_pcc'
            % choose k just higher than the length of small sequence for
            % best precision. (larger k though could increase speed)
            comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))),theoryStruct.isLinearTF,numPixelsAroundBestTheoryMask);
        case 'dtw'
            
            import SignalRegistration.ucr_dtw_score;
            comparisonFun = @(X1,X2, bitX1, bitX2, w) ucr_dtw_score(X1,X2, bitX1,bitX2,w);

%             comparisonStructure{idx}.ucr = ucr_dtw_score(theoryStruct.filename, barC, barB, sets);
        case 'mp'
            import mp.mp_dist_stomp_with_masks;
            comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_dist_stomp_with_masks(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))),[],theoryStruct.isLinearTF);
        case 'mpnan'
            % removes all placements where at least one pixel is taken from a
            % nan region.
            % todo: this function itself does not include possibility of
            % being circular, so inputed X1, X2 should include extra
            % values..
            import mp.mp_profile_stomp_nan_dna;
            comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_profile_stomp_nan_dna(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))),numPixelsAroundBestTheoryMask);

        case 'mpAll' % mp for all subfragments separately
            import mp.mp_dist_stomp_with_masks_all;
             comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_dist_stomp_with_masks_all(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))));

%             comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_dist_stomp_with_masks_all(X1',[X2 X2(1:w-1)]', bitX1', bitX2', w,2^(4+nextpow2(length(X1))));

%             [scoreMatPCC,scorePos,orientation,secondPos, abMatchTemp] = mp_dist_stomp_with_masks(X1,X2, bitX1, bitX2, w, kk,par,islinear)

%             comparisonFun = @(x,y,z,w) unmasked_MP(y,x,z,2^(4+nextpow2(length(x))));
% 
%             error('not yet implemented');
        case 'SEC-C'
            % we here
            error('Check https://github.com/Naderss/SEC_C for implementation');
        case 'awarp'

        case 'hmm'
            % here we stor the hmm method used for structural variation
            % detection. However, the outputs need to be aligned...
            error('not yet implemented..');

        case 'mpdist'
            error('not yet implemented');
        otherwise
            error('undefined method');
    end
            

    % stretch factors
%     stretchFactors = sets.theory.stretchFactors;
%     try % some versions of sets might not have this number, therefore we add a test here
%         minLength = sets.comparison.minLength;
%     catch
%         minLength = 20;
%     end
    
    % number of pixels around best theory match coefficient
%     try
%         numPixelsAroundBestTheoryMask = sets.bitmasking.numPixels;
%     catch
%         numPixelsAroundBestTheoryMask = 0;
%     end

    
    % load theory barcode txt file. For UCR DTW (c++ code), we only need the name of the
    % file so this can be skipped.
    if ~isfield(theoryStruct,'rawBarcode') % enable possibility
        fileID = fopen(theoryStruct.filename,'r');
        formatSpec = '%f';
        theorBar = transpose(fscanf(fileID,formatSpec));
        fclose(fileID);
        try
            fileID = fopen(strrep(theoryStruct.filename,'_barcode','_bitmask'),'r');
            theorBit = transpose(fscanf(fileID,'%f'));
            fclose(fileID);
        catch
            theorBit = ones(1,length(theorBar));
        end
    else
        theorBar = theoryStruct.rawBarcode;
        if ~isempty(theoryStruct.rawBitmask)
            theorBit = theoryStruct.rawBitmask;
        else
            theorBit = ones(1,length(theorBar));
        end
    end
        
    % if we want to take a subtheory
    if nargin >=7
        if ~isempty(positionsOnThry)
            theorBar = theorBar(positionsOnThry(1):positionsOnThry(2));
            theorBit = theorBit(positionsOnThry(1):positionsOnThry(2));
        end
    end

   	import CBT.Hca.Core.filter_barcode; % in case we need to filter barcode

    theorBar = filter_barcode(theorBar,sets);
    
    rezMaxM = cell(1,length(barcodeGen));
    bestBarStretch = zeros(1,length(barcodeGen));
    bestLength = zeros(1,length(barcodeGen));
    % for all the barcodes run
    % parfor
    for idx=1:length(barcodeGen)
%         idx
        % xcorrMax stores the  maximum coefficients
        xcorrMax = zeros(1,length(stretchFactors));
        
        % rezMaz stores the results for one barcode
        rezMax = cell(1,length(stretchFactors));
       
        try
        % barTested barcode to be tested
        barTested = barcodeGen{idx}.rawBarcode;
        catch
        barTested = barcodeGen{idx}.barcode;
        end
        
        barTested = filter_barcode(barTested,sets);

        % in case barcode should be filtered
%         barTested = filter_barcode(barTested, sets.filterSettings);

        % length of this barcode
        lenBarTested = length(barTested);
        
        % barBitmask - bitmask of this barcode
        try
        barBitmask = barcodeGen{idx}.rawBitmask;
        catch
        barBitmask = barcodeGen{idx}.bitmask;
        end

        % run the loop for the stretch factors
        for j=1:length(stretchFactors)
            % here interpolate both barcode and bitmask 
            try
                % simplify scaling: one function (but mask has to be
                % continuous
%                 barC = imresize(barTested(logical(barBitmask)),'Scale' ,[1 stretchFactors(j)]);
                barC = interp1(barTested, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));
                barB = barBitmask(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));

                [rezMax{j}.maxcoef,rezMax{j}.pos,rezMax{j}.or,rezMax{j}.secondPos,rezMax{j}.lengthMatch,~] = comparisonFun(barC, theorBar, barB,theorBit,w);
            catch
                rezMax{j}.maxcoef = zeros(1,3);rezMax{j}.pos=zeros(1,3);rezMax{j}.or=zeros(1,3);rezMax{j}.secondPos=zeros(1,3);rezMax{j}.lengthMatch=0;rezMax{j}.dist=zeros(1,3);
            end
            xcorrMax(j) = rezMax{j}.maxcoef(1);
        end

        % find which stretching parameter had the best score
        [value,b] = max(xcorrMax);

        % select the results for this best stretching parameter and output
        % them. If there were no values computed for this barcode, we don't
        % save anything.
        if ~isnan(value)
            rezMax{b}.bestBarStretch = stretchFactors(b);
            rezMaxM{idx} = rezMax{b};
            bestBarStretch(idx) = stretchFactors(b);
            bestLength(idx) = round(lenBarTested*stretchFactors(b));
%             rezMaxM{idx}.bestBarStretch = bestBarStretch;
        else
            rezMaxM{idx} = rezMax{b};
%              rezMaxM{idx} =nan;
            bestBarStretch(idx) = nan;
            bestLength(idx) = nan;   
        end
    end 
end
% % 
% %% test pccs
% pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
% bar1 = barTested; % or stretched if stretch factor is not 1
% bar2 = theorBar;
% bit1 = barBitmask;
% idxpos = rezMaxM{1}.secondPos(1);
% pos = rezMaxM{1}.pos(1);
% 
% if ~rezMaxM{1}.or
%     frag1 = flipud(bar1);
% else
%     frag1 = bar1;
% end
% % 
% % idxpos = idxpos+1
% frag2 = bar2(pos:pos+length(bar1)-1);
% % 
% bit1 = find(bit1,1,'first');
% % 
% idxpos = idxpos+bit1-1;
% frag1subseq = frag1(idxpos:idxpos+w-1);
% frag2subseq = frag2(idxpos:idxpos+w-1);
% % 
% pcc(frag1subseq',frag2subseq')
