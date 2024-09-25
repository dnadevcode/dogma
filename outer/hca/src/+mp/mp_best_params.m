function [ maxcoef,pA,or,pB ] = mp_best_params( mp, mpI, mpD, numBestPar, mask)
    % get_best_parameters_mp
    %
    % todo: merge to get_best_parameters, so that all cases are in the same
    %
    %   Args:
    %
    %   Returns:
    %
    % input xcorrs - cross correlation coefficients
    % numBestPar - number of coefficients to output
    % output rezMax
    
    if nargin < 4
        numBestPar = 3; % default is 3
    end
    
    if nargin < 5
        mask = 0;
    end
     
    % f gives the max over the two rows, s gives the indice of
    % which row has the max (1 or 2). TODO: maybe change this to be a match
    % table instead - corresponds to structural variations better.
    maxcoef = zeros(1,numBestPar);
    or = zeros(1,numBestPar);
    pA = zeros(1,numBestPar);
    pB = zeros(1,numBestPar);

    for ii=1:numBestPar
        [f, s] = max(mp);

        % sort the max scores, ix stores the original indices
        [ b, ix ] = sort( f(:), 'descend','MissingPlacement','last' );

        % choose the best score
        indx = b(1) ;
        % save the best max score and orientation
        maxcoef(ii) = indx;
        
        pA(ii) = s(ix(1));
    
        % now orientation is easy!
        or(ii) = mpD(s(ix(1)))+1;
        
        pB(ii) = mpI(s(ix(1)));

        % finally, save the position. This can have two cases,
        % depending on the value of s
        mp(s) = nan;
        
        mp(logical((mpI > mpI(s(ix(1)))-mask).*(mpI < mpI(s(ix(1)))+mask))) = nan;

    end

        
end