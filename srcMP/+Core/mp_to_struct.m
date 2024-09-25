function [overlapStruct] = mp_to_struct(mp1,mpI1,baridx2,stridx,w, sF,barStruct)
    
    % Args: mp1 - matrix profile, mpI1 - matrix profile index, baridx2 -
    % barcode indexes, stridx - stretch factors and orientations, w -
    % overlap window width, barStruct - barcode structure

    % Returns:  overlapStruct - overlap structure
    
    N = length(mp1);
    
    import Core.get_full_overlap_score;

    % pre-define all fields
    overlapStruct = struct('score', cell(N,N), 'fullscore', cell(N,N), 'overlaplen', cell(N,N), 'pA', cell(N,N),...
        'pB', cell(N,N), 'or', cell(N,N), 'bestBarStretch', cell(N,N),'lenA', cell(N,N),...
        'lenB', cell(N,N));
    %     overlapStruct = [];%cell(1,length(mp1));
    
    for k=1:length(mp1) % loop through mp
        for iy=1:length(mp1)
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
            overlapStruct(k,iy).partialScore  = nan;
            overlapStruct(k,iy).partialLength  = nan;

        end
    end


for k=1:length(mp1) % loop through mp

    posfirstAll = zeros(2,length(sF)); % first pos in re-scale factor. could be already given when generating data.
    for ii=1:length(sF)
        posfirstAll(1,ii) = find(stridx{k}==ii,1,'first'); % if we take just the best
        posfirstAll(2,ii) = find(stridx{k}==-ii,1,'first'); % if we take just the best
    end


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
            overlapStruct(k,iy).h = w;
            % this can be done outside to speed things up, we might want
            % this only for some of the barcodes
            [overlapStruct(k,iy).fullscore,overlapStruct(k,iy).overlaplen, overlapStruct(k,iy).lenB , overlapStruct(k,iy).lenA,overlapStruct(k,iy).partialScore,...
                 overlapStruct(k,iy).partialLength] = get_full_overlap_score(overlapStruct(k,iy).pA,overlapStruct(k,iy).pB,...
                 overlapStruct(k,iy).bestBarStretch, overlapStruct(k,iy).or,barStruct([k iy]),w);

        end
    end
end

end