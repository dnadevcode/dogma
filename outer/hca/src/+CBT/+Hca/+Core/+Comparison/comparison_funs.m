function [comparisonFun] = comparison_funs(comparisonMethod,isLinearTF,numPixelsAroundBestTheoryMask)


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
    
end

