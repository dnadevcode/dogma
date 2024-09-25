function [ rezMaxM,bestBarStretch,bestLength ] = on_compare_rescale_barcode(barcodeGen,theorBar, theorBit, isLinearTF, comparisonMethod,stretchFactors,w,numPixelsAroundBestTheoryMask,sigma1)
    % on_compare_barcode / in this theorBar opening is moved to the outside
    % loop to speed things up possibly
    % Compares experiments to single theory
    %     Args:
    %         sets: settings structure
    %         barcodeGen: barcode structure
    %         theorBar - theory barcode,theorBit - theory bitmask
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
            comparisonFun = @(x,y,z,w,u) unmasked_MASS_DOT_CC(y,x,z,w,2^(4+nextpow2(length(x))), isLinearTF,numPixelsAroundBestTheoryMask);            
       case 'masked_pcc'
            % choose k just higher than the length of small sequence for
            % best precision. (larger k though could increase speed)
            comparisonFun = @(x,y,z,w,u) masked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))), isLinearTF,numPixelsAroundBestTheoryMask);

        case 'mass_pcc'
            % choose k just higher than the length of small sequence for
            % best precision. (larger k though could increase speed)
            comparisonFun = @(x,y,z,w,u) unmasked_MASS_PCC(y,x,z,w,2^(4+nextpow2(length(x))), isLinearTF,numPixelsAroundBestTheoryMask);
        case 'dtw'
            
            import SignalRegistration.ucr_dtw_score;
            comparisonFun = @(X1,X2, bitX1, bitX2, w) ucr_dtw_score(X1,X2, bitX1,bitX2,w);

%             comparisonStructure{idx}.ucr = ucr_dtw_score(theoryStruct.filename, barC, barB, sets);
        case 'mp'
            import mp.mp_dist_stomp_with_masks;
            comparisonFun = @(X1,X2, bitX1, bitX2, w, kk) mp_dist_stomp_with_masks(X1',X2', bitX1', bitX2', w,2^(4+nextpow2(length(X1))),[],isLinearTF);
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

        case 'hmm'
            % here we stor the hmm method used for structural variation
            % detection. However, the outputs need to be aligned...
            error('not yet implemented..');

        case 'mpdist'
            error('not yet implemented');
        otherwise
            error('undefined method');
    end
              

   	import CBT.Hca.Core.filter_barcode; % in case we need to filter barcode

    rezMaxM = cell(1,length(barcodeGen));
    bestBarStretch = zeros(1,length(barcodeGen));
    bestLength = zeros(1,length(barcodeGen));
    % for all the barcodes run
    % parfor
    for idx=1:length(barcodeGen)
        idx
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
            % here we perform the rescale:
            
            len1 = round(lenBarTested*stretchFactors(j));
            intPts = linspace(1,len1,lenBarTested);
            barC = interp1(intPts, barTested,1:len1,'spline','extrap');
            barB = barBitmask(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));
            
            % need sigma
            sigma1 = 300/169.2; % later add vias ettings
            sigma2=sigma1*len1/lenBarTested; % it should differ by lenBar1/lenBarTested

            % how much is sigma rescaled by interpolation. Could also check how much
            % higher order approximation would improve things. How much the accuracy
            % depends on sigma?
            if sigma2 > sigma1
                barC = run_deconvolve(barC,sqrt(sigma2.^2-sigma1.^2));
            else
                if sigma1 > sigma2
                    barC = run_convolve(barC,sqrt(sigma1.^2-sigma2.^2));
                end
            end
    
%             barC = interp1(barTested, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));
            try
                [rezMax{j}.maxcoef,rezMax{j}.pos,rezMax{j}.or,rezMax{j}.secondPos,rezMax{j}.lengthMatch,~] = comparisonFun(barC, theorBar, barB,theorBit,w);
            catch
                rezMax{j}.maxcoef = 0;rezMax{j}.pos=0;rezMax{j}.or=0;rezMax{j}.secondPos=0;rezMax{j}.lengthMatch=0;rezMax{j}.dist=0;
            end
            xcorrMax(j) = rezMax{j}.maxcoef(1);
        end

        % find which stretching parameter had the best score
        [value,b] = max(xcorrMax);

        % select the results for this best stretching parameter and output
        % them. If there were no values computed for this barcode, we don't
        % save anything.
        if ~isnan(value)
            rezMaxM{idx} = rezMax{b};
            rezMaxM{idx}.rescaleFactor = stretchFactors(b);
            bestBarStretch(idx) = stretchFactors(b);
            bestLength(idx) = round(lenBarTested*stretchFactors(b));
            rezMaxM{idx}.allMax = cellfun(@(x) x.maxcoef(1),rezMax);
            rezMaxM{idx}.allPos = cellfun(@(x) x.pos(1),rezMax);

%             rezMaxM{idx}.bestBarStretch = bestBarStretch;
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
