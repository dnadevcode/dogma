%     bTS = rand(1000,1);
%     aTS = rand(100,1);
%     r = 50;
%     kk =128;
% [mp,mpI] = mp_profile(aTS, bTS,r,kk);

function [maxcoef, pB, or, pA, r, mp, mpI] = mp_profile_stomp_nan_dna(aTS, bTS, bitA, maskB, r, kk,numPixelsAroundBestTheoryMask)

    if nargin < 7
        numPixelsAroundBestTheoryMask = 50;
    end
    % mp profile stomp for DNA comparison (forward and reverse strand)
    % bitmasks should already included in aTS and bTS as NANs
    mpLength = length(aTS);%-r+1;
    
%     mp2 = zeros(mpLength,1); % instead of making twice the size, report a 
%     mpI2 = zeros(mpLength,1); % a negative number in mpI? Then all procedure simplifies.
    
    % mask both aTS and bTS (bTS is masked if it includes at least 50
    % undefined base-pairs
%     nanValues = movmean(maskB,length(aTS),'Endpoints','discard');
    bTS(maskB > 50) = nan;        % 50 is estimate.., mean approx 50 bp  (1/10th) in a pixel were nan's
    aTS(~bitA) = nan;
    
    % forward and backward, together
    [mp, mpI] = mp_profile_stomp([aTS; nan; aTS(end:-1:1)], bTS, r, kk);

    [maxcoef,pos] = max(mp);
    pB = mpI(pos);

    pA = mod(pos, mpLength+1);
    or = (pos > mpLength)+1; % now the output is simple. Shows position on barB (pB), barA (pA) and orientation of barA


    % backward, simple

    % reverse, compare reversed aTS vs bTS, also save in reversed,
    % so mp2(1) maps flip[aTS(1:r)] to bTS, and we compare this to aTS(1:r)
    % vs bTS
    %[mp2(end:-1:1), mpI2(end:-1:1)] = mp_profile_stomp(aTS(end:-1:1), bTS,r,kk);
    
    % update mp positions where reverse has higher PCC than direct
	% mpPos(:,1) = mp(:) < mp2;
    % mp(mpPos) = mp2(mpPos);
    % mpI(mpPos) = mpI2(mpPos); % don't need to add anything, since we save updatePos also
    
    % mpI is now position of the subfragment. But instead we might want to
    % compare the actual position // maybe not all, just best?
    %     import mp.convert_subfragment_pos_to_query;
    %     [mpPIQ,mpI] = convert_subfragment_pos_to_query(mpI,length(aTS),r,mpPos);
    
    % now get parameters for best positions: maxcoef, pos, or, idxpos
    % (which means position on second barcode)
    % consider how this would be extracted from SCAMP also.
% 	import mp.mp_best_params; % idxpos should be best pos on second barcode
%     [ maxcoef, pB, or, pA ] = mp_best_params( mp, mpI, mpPos, 5, numPixelsAroundBestTheoryMask);    
end

function [mp,mpI] = mp_profile_stomp(aTS, bTS, r, kk)
    % simple method to compute matrix profile for time series data comping
    % from DNA barcoding
    %
    % Args:
    %   aTS, bTS, r
    %
    % Returns:
    %   mp: matrix profile,
    %   mpI: matrix profile index
    
    mpLength = length(aTS)-r+1;
    mpColLength = length(bTS)-r+1;

    % want to define nan on all the values affected by nanMaskbTS, this
    % will be nanMaskbTS+relements before
    nansB = isnan(bTS); % mask bTS
    bTS(nansB) = 1;
    nansA = isnan(aTS); % mask aTS
    aTS(nansA) =1;
  
    [X, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(bTS, r);
    [Y, sumy2, sumy, meany, sigmay2, sigmay] = fastfindNNPre(aTS, r);
    
    % check if this is correct // nan should be for all that 
    meany(nansA(1:length(meany))) = nan;
    sigmay(nansA(1:length(meany))) = nan;
%     firstNanA = find(nansA,1,'first');
    
    % find all 0-1 sequences.
    nanStarts = strfind(nansA',[0 1]);
%     nanStarts = nanStarts(nanStarts<=length(aTS)-r); % is it n-elements
%     forward or backward that have to be bitmasked?
    for i=1:length(nanStarts)
        meany(max(1,nanStarts(i)-r+2):min(nanStarts(i),length(meany))) = nan;
        sigmay(max(1,nanStarts(i)-r+2):min(nanStarts(i),length(meany))) = nan;
    end

    %      find(nansA,[0;1])
%     for i=1:length(nanStarts)
%         meany(nanStarts(i)-r+1:nanStarts(i)-1) = nan;
%         sigmay(nanStarts(i)-r+1:nanStarts(i)-1) = nan;
%     end

    
%     firstNan = find(nansB,1,'first');
    nanStarts = strfind(nansB',[0 1]);
%     nanStarts = nanStarts(nanStarts<=length(bTS)-r);

    meanx(nansB(1:length(meanx))) = nan;
    sigmax(nansB(1:length(meanx))) = nan;
    
    for i=1:length(nanStarts)
        meanx(max(1,nanStarts(i)-r+2):min(nanStarts(i),length(meanx))) = nan;
        sigmax(max(1,nanStarts(i)-r+2):min(nanStarts(i),length(meanx))) = nan;
    end

    qRow = zeros(mpLength,1); % same as MP
    qCol = zeros(mpColLength,1); % columns distance profile

    % Start from
    dist = FAST_CC(bTS,aTS(1:r), kk );
    dist2 = FAST_CC(aTS,bTS(1:r), kk );

    qCol(:)= ((dist./r)-meanx*meany(1))./(sigmax.*sigmay(1));
    qRow(:) = ((dist2./r)-meanx(1)*meany)./(sigmax(1).*sigmay);
    
    mp = zeros(mpLength,1); % instead of making twice the size, report a 
    mpI = zeros(mpLength,1); % a negative number in mpI? Then all procedure simplifies.
    updatePos=false(mpLength,1);

    % evaluate initial matrix profile
    mp(:) = qRow;
    mpI(:) = 1;
    
    updatePos(:,1) = mp(:) < qRow;
    mp(updatePos) = qRow(updatePos);
    mpI(updatePos) = mpColLength+1;
    [ mp(1), mpI(1)] = max(qCol);

    for i = 2:mpLength
        dist(2:mpColLength) = dist(1:mpColLength-1)- bTS(1:mpColLength-1).*aTS(i-1)+bTS(r+1:end).*aTS(i+r-1);

        % update first values
        dist(1) = dist2(i);
        qCol= ((dist./r)-meanx*meany(i))./(sigmax.*sigmay(i));
        [ mp(i), mpI(i)] = max(qCol);
    end
end


function [dist] = FAST_CC(x, y, k)
    % FAST CC- batch processing of CC on short segments 
    %
    %   Args:
    %       x, y, k
    %
    %   Returns:
    %       dist, distance matrix
    
    %x is the data, y is the query
    m = length(y);
    n = length(x);
    mpLen = n-m+1;
    dist = zeros(mpLen,1);

    y = y(end:-1:1); %Reverse the query
    y(m+1:k) = 0; %append zeros

    Y = fft(y);

    for j = 1:abs(k-m+1):n-k+1

        %The main trick of getting dot products in O(n log n) time
        X = fft(x(j:j+k-1));

        Z = X.*Y;
        z = ifft(Z);
        dist(j:j+k-m) = z(m:k);
    end

     if isempty(j)
        j = 0; % if nothing was computed
        k = n;
     else
        j = j+k-m;
        k = n-j; % number of points left
     end
     
    if k >= m % if k < m, there are not enough points on long barcode to compute more PCC's
        
        %The main trick of getting dot products in O(n log n) time
        X = fft(x(j+1:n));

        y(k+1:end)= [];

        Y = fft(y);
        
        Z = X.*Y;
        z = ifft(Z);

        dist(j+1:n-m+1) = z(m:k);
    end    
end


% m is winSize
function [X, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
    n = length(x);
    x(n+1:2*n) = 0;
    X = fft(x);
    cum_sumx = cumsum(x);
    cum_sumx2 =  cumsum(x.^2);
    sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
    sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
    meanx = sumx./m;
    sigmax2 = (sumx2./m)-(meanx.^2);
    sigmax = sqrt(sigmax2);
end
