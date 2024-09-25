function [overlapStruct] = mp_to_struct_multi(N, mp1, mpI1, baridx2, stridx, w, sF)
    
    % Args: mp1 - matrix profile, mpI1 - matrix profile index, baridx2 -
    % barcode indexes, stridx - stretch factors and orientations, w -
    % overlap window width, barStruct - barcode structure

    % Returns:  overlapStruct - overlap structure.
    
%     N = length(mp1);
    
    import Core.get_full_overlap_score;

    % pre-define all fields
    overlapStruct = struct('score', cell(1,N), 'fullscore', cell(1,N), 'overlaplen', cell(1,N), 'pA', cell(1,N),...
        'pB', cell(1,N), 'or', cell(1,N), 'bestBarStretch', cell(1,N),'lenA', cell(1,N),...
        'lenB', cell(1,N));


    % position after length re-scalings
    posfirstAll = zeros(2,length(sF)); % first pos in re-scale factor. could be already given when generating data.
    for ii=1:length(sF)
        posfirstAll(1,ii) = find(stridx == ii,1,'first'); % if we take just the best
        posfirstAll(2,ii) = find(stridx == -ii,1,'first'); % if we take just the best
    end

    
    for k=1:N % loop through mp
        posA = baridx2{1}==k;
        subMP =  mp1{1}(posA);
        subMPI = mpI1{1}(posA);
        [mpS,locs] = sort(subMP,'descend','MissingPlacement','last');
        overlapStruct(1,k).score = mpS(1);
        overlapStruct(1,k).pB = locs(1); % position on root (un-rescaled barcode)

        % orientation
        srtIdx =  stridx(subMPI(locs(1))+1);
        overlapStruct(1,k).or = -0.5*sign(srtIdx)+3/2;

        posfirst = posfirstAll(  overlapStruct(1,k).or, abs(srtIdx)); %find(stridx{k}==srtIdx,1,'first'); % if we take just the best

        overlapStruct(1,k).pA = subMPI(locs(1))+1-posfirst+1;
        overlapStruct(1,k).bestBarStretch = sF(abs(srtIdx));
        overlapStruct(1,k).h = w;

        % this can be done outside to speed things up, we might want
        % this only for some of the barcodes. For circular overlaps not
        % exactly correct
%         [overlapStruct(k,iy).fullscore,overlapStruct(k,iy).overlaplen, overlapStruct(k,iy).lenB , overlapStruct(k,iy).lenA,overlapStruct(k,iy).partialScore,...
%              overlapStruct(k,iy).partialLength] = get_full_overlap_score(overlapStruct(k,iy).pA,overlapStruct(k,iy).pB,...
%              overlapStruct(k,iy).bestBarStretch, overlapStruct(k,iy).or,barStruct([k iy]),w);

    end
end