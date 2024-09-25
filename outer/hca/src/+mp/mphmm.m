%     bTS = rand(1000,1);
%     aTS = rand(100,1);
%     r = 50;
%     kk =128;
% [mp,mpI] = mp_profile(aTS, bTS,r,kk);

function [maxcoef,pos,or,idxpos, mp,mpI,mpPos,mpPIQ] = mphmm(aTS, bTS, r, kk)
    % mp profile stomp for DNA comparison (forward and reverse strand)
    % bitmasks should already included in aTS and bTS as NANs
    mpLength = length(aTS)-r+1;
    
    mp2 = zeros(mpLength,1); % instead of making twice the size, report a 
    mpI2 = zeros(mpLength,1); % a negative number in mpI? Then all procedure simplifies.
    
    % if aTS assumed linear, here we create
    % start simple case: no nan
%     aTS = [aTS; nan; flipud(aTS)];
    % if circular:
    aTS = [aTS;aTS(1:r-1)];
    
%     figure,plot(aTS)
    
    % todo: merge forward and reverse into single run, adding a "nan" in
    % between, so matches combining both would be forbidden
    % forward, simple
    [mp,mpI, path] = mp_profile_stomp(aTS, bTS, r, kk);
    
    
    % not sure about the following, still need to fix for circular/
    % reversed data..
    % can we convert path to table matrix? 
    
    %% NO need to compute reverse as this is already included in aTS
    % reverse, compare reversed aTS vs bTS, also save in reversed,
    % so mp2(1) maps flip[aTS(1:r)] to bTS, and we compare this to aTS(1:r)
    % vs bTS
%     [mp2(end:-1:1),mpI2(end:-1:1)] = mp_profile_stomp(aTS(end:-1:1), bTS,r,kk);
    
    % update mp positions where reverse has higher PCC than direct
	mpPos(:,1) = mp(:) < -inf;
%     mp(mpPos) = mp2(mpPos);
%     mpI(mpPos) = mpI2(mpPos); % don't need to add anything, since we save updatePos also
    
    % mpI is now position of the subfragment. But instead we might want to
    % compare the actual position // maybe not all, just best?
    import mp.convert_subfragment_pos_to_query;
    [mpPIQ] = convert_subfragment_pos_to_query(mpI,length(aTS),r,mpPos);
    
    % now get parameters for best positions: maxcoef, pos, or, idxpos
    % (which means position on second barcode)
    % consider how this would be extracted from SCAMP also.
	import mp.mp_best_params; % idxpos should be best pos on second barcode
    [ maxcoef,pos,or,idxpos ] = mp_best_params( mp, mpI, mpPIQ, mpPos, 3, 50);    
end

function [mp,mpI,path] = mp_profile_stomp(aTS, bTS, r, kk)
    % simple method to compute matrix profile for time series data comping
    % from DNA barcoding
    %
    % Args:
    %   aTS, bTS, r
    %
    % Returns:
    %   mp: matrix profile,
    %   mpI: matrix profile index
    %   path: computed Viterbi path
    
    %TODO: versions for both AB and BA MP/ maybe one of them can be run more
    %efficiently
    
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
    
%     %% NAN STUFF
%     meany(nansA(1:length(meany))) = nan;
%     sigmay(nansA(1:length(meany))) = nan;
%     
%     % find all 0-1 sequences.
%     nanStarts = strfind(nansA',[0 1]);
%     nanStarts = nanStarts(nanStarts<=length(aTS)-r);
%     for i=1:length(nanStarts)
%         meany(max(1,nanStarts(i)-r+2):nanStarts(i)) = nan;
%         sigmay(max(1,nanStarts(i)-r+2):nanStarts(i)) = nan;
%     end
% 
%     nanStarts = strfind(nansB',[0 1]);
%     nanStarts = nanStarts(nanStarts<=length(bTS)-r);
% 
%     meanx(nansB(1:length(meanx))) = nan;
%     sigmax(nansB(1:length(meanx))) = nan;
%     
%     for i=1:length(nanStarts)
%         meanx(max(1,nanStarts(i)-r+2):nanStarts(i)) = nan;
%         sigmax(max(1,nanStarts(i)-r+2):nanStarts(i)) = nan;
%     end
    
    % end nan stuff

    qRow = zeros(mpLength,1); % same as MP
    qCol = zeros(mpColLength,1); % columns distance profile

    % Start from
    dist = FAST_CC(bTS,aTS(1:r), kk );
    dist2 = FAST_CC(aTS,bTS(1:r), kk );

    % calculate row and column
    qCol(:)= ((dist./r)-meanx*meany(1))./(sigmax.*sigmay(1));
    qRow(:) = ((dist2./r)-meanx(1)*meany)./(sigmax(1).*sigmay);
    
    % maybe plot qCol and qRow for write-up example
    
    % output: we compute both mp and mpI, and the hmm. Time consumption is
    % similar...
    
    mp = zeros(mpLength,1); % instead of making twice the size, report a 
    mpI = zeros(mpLength,1); % a negative number in mpI? Then all procedure simplifies.
    path = [];
    updatePos=false(mpLength,1);

    % evaluate initial matrix profile
    mp(:) = qRow;
    mpI(:) = 1;
    
    updatePos(:,1) = mp(:) < qRow;
    mp(updatePos) = qRow(updatePos);
    mpI(updatePos) = mpColLength+1;
    [ mp(1), mpI(1)] = max(qCol);
    
    %% HMM
    enableViterbi = 1;
    if enableViterbi
        % qCol could be converted to p-values here - for the use in HMM
        % INITIALIZE HMM
        gapScore = 0.85; % we add a score for having a gap. We can move to a gap state, and we can 
        % define intial scores
        initial = [qCol; nan; gapScore];
        [delta, psi, pM] = setup_data(initial, mpLength);
        lastElt = ((aTS(r)-meany(1))./sigmay(1) -(bTS(r:end)-meanx)./sigmax).^2;
        delta = delta';
        psi = psi';
    end
%     f = figure,plot(delta(:,1))
%     saveas(f,'ex3.eps','epsc');
%     
    harmonicMean = @(x,y,j) j./((j-1)./x+1./y); 

%     % should work similar also in case of allowing circular data.
% we could compute both row and column at the same time, going from top
% left to bottom right
    for i = 2:mpLength
        dist(2:mpColLength) = dist(1:mpColLength-1)- bTS(1:mpColLength-1).*aTS(i-1)+bTS(r+1:end).*aTS(i+r-1);

        % update first values
        dist(1) = dist2(i);
        qCol = ((dist./r)-meanx*meany(i))./(sigmax.*sigmay(i));
        
        % als
%         ed2Col = 2*r*(1-qCol);
        
        %% Viterbi
        if enableViterbi
            % In this version, those scores that are less than gapScore, will
            % move to gap state. 
%             prod
            % last element score addition, this somehow contributes the
            % moving probability, and can be combined together with MP
            % score
            
            % we can choose to add this last element contribution to the
            % algorithm. Also can do this for last two elements, or last
            % three elements, or last n-elements. If this score is too
            % high, we can discard safely... useful??
            
            % instead of expanding and using lower bound, we compute the
            % last element of each sequence / numeric issues?
            lastElt = ((aTS(i+r-1)-meany(i))./sigmay(i) -(bTS(r:end)-meanx)./sigmax).^2;
%             lastElt = (bTS(r+1:end).*aTS(i+r-1)./r-meanx*meany(i))./(sigmax.*sigmay(i));
%            (( aTS(51)-meany(2))/sigmay(2) - (bTS(101)-meanx(51))/sigmax(51)).^2

            % here we have both dist score has to be smaller than
            % something, and also last element score. Note that the last
            % element score is a lower bound on  the total distance
            % most important part: deciding when the current path goes to
            % gap state
            pass = find(qCol > gapScore); % gap score depend on what kind of random data we have?
            % hyper-parameter
            % could be that all are nan's, then inf should be larger
            [ delta(end-1,i), idx ] = max([-inf; qCol(pass)]);

            % add previous score, which is taken from one element on the left,
            % and only if idx>1, i.e. there are elements passing threshold
            if idx > 1
%                 delta(end-1,i) = delta(end-1,i) + delta(pass(idx-1),i-1);
                delta(end-1,i) = harmonicMean(delta(pass(idx-1),i-1),gapScore,i); % what score should we give/ should not be the same as for the last element
                psi(end-1,i) = pass(idx-1)-1;
                
                % delta(end-1,i) - score for subsequence Q_{i,r} vs D_{j,r}
                % where j = pass(idx-1)
                % now we want to combine 
                % $dist(Q_{i,r}, D_{j,r})$ and the distance
                % $dist(Q_{i-1,r}, D_{j-1,r})$
                % 
                % we know that  delta(pass(idx-1),i-1), delta(end-1,i) - 
                % here  delta(pass(idx-1),i-1) is combination of previous
                % i-1 scores.
                
                % One way: just use harmonic mean
                % we have harmonic mean of i-1 elements:
                % i/sum_k(1/dist_k)
                % now, we add dist_k+1
                % so we want
                % (i+1)/sum_k(1/dist_k)
            end

            % now other states/ just add the previous probability // needs an
            % example %todo: what if one of the elements is INF?
            delta(2:end-2,i) = harmonicMean(delta(1:end-3,i-1),qCol(2:end),i);
            psi(2:end-2,i) =  psi(1:end-3,i-1); % and these are also one before

            % now gap. 
            [delta(end,i), gapSt] = max([harmonicMean(delta(end,i-1),gapScore,i), delta(end-1,i)]);

            if gapSt == 2
                psi(end,i) = mpColLength+1;
            end
        end
        % end MP stuff
        
        [ mp(i), mpI(i)] = max(qCol);
    end
    
    
    
    % Viterbi traceback
    if enableViterbi
        [score, idx]  = max(delta(:,end));

        % traceback..
        path = [mpLength, idx];

        for k= mpLength:-1:2
            newV = psi(idx,k);
            path = [path; k, newV];
            if newV == mpColLength+1
                newV = psi(newV,k);
                path = [path; k, newV];
            end
            newV = idx;        
        end
    end
    
    % forward/backward?
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




    function [delta, psi,pM] = setup_data(scores, L, rM, rG)
        % set up HMM-MP data
    
        if nargin < 3
%             rM = 16.7;
%             rG = 1;
%             rM = 0.5;
%             rG = 0.5;

        end
        
%         p_MM = rM/(1+rM); % If we allowed movement to more states, this would be modified
%         p_MG = 1/(1+rM);
%         p_GM = 1/(rG + length(scores)-2);
%         p_GG = rG/(rG +length(scores)-2);
        
        
        % probabilities (in log10) of match to match, match to gap, gap to
        % match, gap to gap. 
%         pM = log10([ p_MM p_MG; p_GM p_GG]);
%         pM = [ p_MM p_MG; p_GM p_GG];
        % remove transition scores for hmm-mp..
        pM = zeros(2,2);

        % Alternatively we compare p-value scores (PCCs) p1 and p2.
        % might be that we include a lot of pixels at the edges that we don't want (
        % especially for long initial length)

  
        % emission of gap. this is closely related to score when comparing
        % non-match fragments or partially overlapping fragments)
        emG = log10(0.5); 

   
        N = length(scores); % total number of states

%         p = log10(1/(N-1)); % each state has prob 1/(2M+1)

        % delta, psi, and xi, which we keep track of throughout the algorithm.
        % how to reduce large memory use here? maybe just keep track of psi
        % matrix rather than delta?
        delta = -inf(L,N);
        psi = zeros(L,N);

        % initialization (if M is large, here we might have memory issues for
        % storing a very big matrix)
        % so [S columns forward, S columns backward, and Gap]

        % choose to use either col/row operations. Could use more faithful
        % sumation to reduce rounding errors
        delta(1,:) = scores;
%         delta(1,:) = [scorefun(data(1),query)+p scorefun(data(1),flipud(query))+p -inf p+emG]; % so for reverse, the states go M,M-1,M-2,etc
        psi(1,:) = [1:N]; % initialize store which match state each pixel corresponds to
        psi(:,N) = N;
%         xi(:) = 1; % all xi states initialized to 1

    end

