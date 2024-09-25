function [overlapStruct] = mp_res_to_struct_mini(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct)

% mini resultsto struct. Allow overlap to be circular
% totSF = length(sF);


overlapStruct = [];%cell(1,length(mp1));
for k=1:length(mp1)

    posfirstAll = zeros(2,length(sF)); % first pos in re-scale factor. could be already given when generating data.
    for ii=1:length(sF)
        posfirstAll(1,ii) = find(stridx{k}==ii,1,'first'); % if we take just the best
        posfirstAll(2,ii) = find(stridx{k}==-ii,1,'first'); % if we take just the best
    end
    origL = length(barStruct(k).rawBarcode(barStruct(k).rawBitmask));


    for iy=1:length(mp1)
        if k ~= iy
            subMP = mp1{k}(baridx2{k}==iy);
            subMPI = mpI1{k}(baridx2{k}==iy);
            [mpS,locs] = sort(subMP,'descend','MissingPlacement','last');
            overlapStruct(k,iy).score = mpS(1);
            overlapStruct(k,iy).pB = locs(1); % position on root (un-rescaled barcode)

            %
            srtIdx =  stridx{k}(subMPI(locs(1))+1);
            overlapStruct(k,iy).or = sign(srtIdx);

            posfirst = posfirstAll(-0.5*overlapStruct(k,iy).or+1.5, abs(srtIdx)); %find(stridx{k}==srtIdx,1,'first'); % if we take just the best

            overlapStruct(k,iy).pA = subMPI(locs(1))+1-posfirst+1;
            overlapStruct(k,iy).bestBarStretch = sF(abs(srtIdx));
            overlapStruct(k,iy).h = MIN_OVERLAP_PIXELS;
            overlapStruct(k,iy).allposA = subMPI(locs)+1 ;
            overlapStruct(k,iy).allposB = locs;


            pA = overlapStruct(k,iy).pA ;
            pB = overlapStruct(k,iy).pB ;
            h = overlapStruct(k,iy).h;
            orr = overlapStruct(k,iy).or ;
            bBar = imresize(barStruct(k).rawBarcode(barStruct(k).rawBitmask),'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
            aBar = barStruct(iy).rawBarcode(logical(barStruct(iy).rawBitmask));


            if orr~=1
                bBar = fliplr(bBar);
            end
%             zscore(aBar(pB:pB+h-1),1)*zscore(bBar(pA:pA+h-1),1)'/h

            % do full overlap
            lpA = length(bBar); lpB = length(aBar);

            st = min(pA,pB); % which barcode is more to the left
            stop = min(lpA-pA+1,lpB-pB+1);

            aFul = aBar(pB-st+1:pB+stop-1);
            bFul = bBar(pA-st+1:pA+stop-1);
            overlapStruct(k,iy).fullscore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
            overlapStruct(k,iy).overlaplen = length(aFul);
            overlapStruct(k,iy).lenB = length(aBar);
            overlapStruct(k,iy).lenA = length(bBar);

            % placements and scores for overlap
%             for locations on mp
            overlapStruct(k,iy).fulloverlapScoreB = subMP([pB-st+1:pB+stop-MIN_OVERLAP_PIXELS]);
            overlapStruct(k,iy).fulloverlapPosRoot = [pB-st+1:pB+stop-MIN_OVERLAP_PIXELS];
            overlapStruct(k,iy).localSFoverlap = stridx{k}(subMPI(overlapStruct(k,iy).fulloverlapPosRoot)+1);

            overlapStruct(k,iy).fulloverlapor = sign(overlapStruct(k,iy).localSFoverlap);

            pos = arrayfun(@(x,y) posfirstAll(-0.5*x+1.5,abs(y)),overlapStruct(k,iy).fulloverlapor,overlapStruct(k,iy).localSFoverlap);
            overlapStruct(k,iy).fulloverlapPosA = subMPI(overlapStruct(k,iy).fulloverlapPosRoot)'+1-pos+1;

            
            overlapStruct(k,iy).fulloverlapPosARescaled = overlapStruct(k,iy).fulloverlapPosA./sF(abs(overlapStruct(k,iy).localSFoverlap));

%             overlapStruct(k,iy).h = MIN_OVERLAP_PIXELS; 
            % plot overlap:
%             import Core.plot_match_simple;
%             [f] = plot_match_simple(barStruct, overlapStruct,k,iy);

        else
            overlapStruct(k,iy).score  = nan;
            overlapStruct(k,iy).fullscore  = nan;
            overlapStruct(k,iy).overlaplen  = nan;
            overlapStruct(k,iy).pA  = nan;
            overlapStruct(k,iy).or  = nan;
            overlapStruct(k,iy).pB  = nan;
            overlapStruct(k,iy).score  = nan;
            overlapStruct(k,iy).lenA  = nan;
            overlapStruct(k,iy).lenB  = nan;
            overlapStruct(k,iy).bestBarStretch  = nan;

        end
    end
end

end