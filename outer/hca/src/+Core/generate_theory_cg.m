function [thry, bitmask] = generate_theory_cg(Y, ts, cutPointsL, cutPointsR, extraL, extraR, ...
     k, m, circular, numPx, sets)

    import CBT.Hca.Core.Theory.compute_theory_wrapper;

    numLoops = length(cutPointsL);
%% THINGS THAT ARE IMPORTANT
% - when applying PSF(Gaussian), convolution for which points is given as an
% output
% - which elements to use for each pixel
% - how to convert each pixel bp value to px value
% - how to deal with circular data
% - how to deal with undefined data
% - what method to use to convert bp to photon counts (i.e. transfer
% matrix/GC content, patter)
% 
% % y has 1:m centered Gaussian. We can show here how it is 
% x = zeros(length(Y),1);
% x(1)=1; x(5000)=2;
% z= ifft(fft(x).*Y);
% % figure,plot(Z) / m/2-m/2+1 is where the Gaussian peak is
% [a,b] = findpeaks(z,'SortStr','descend','NPeaks',2);
% so x(1)=1 now end's up being m/2

% gives convolved values from m/2 to k-1, those to the left are not
% interesting since convolve values circularly.
% So if there was some elements (m/2)+pxSize/2 added in the set up, z(m:k-1) will
% have a correct starting place. We shift it a little bit based on the
% cutPoints place

% data is divided 1:k, k-m+1:2*k-m, 2*k-2*m+1:3*k-2*m, 3*k-3*m+1:4*k-3*m.
% Now, for each one of these, we calculated cutPointsL, cutPointsR

% idx = ;
% this could be changed slightly, if i.e. we allow overlap
% cutPointsL = floor((1:numPx)*pxSize-pxSize+1);
% cutPointsR = floor((1:numPx)*pxSize);

% floor((i*pxSize)+pxSize); % so we take [-pxSize:pxSize] around the value of interest

% At j'th step, the data we take is
% j*(k-m)+1:k+j*(k-m), 

% j=0, 1:k, but z(m:k-1) 1
%  k-m+1:2*k-m
%  
%  j=5:

% added: m/2+round(pxSize)+1
% count how many loop's there are going to be
% 1:k...  j*(k-m)+1:k+j*(k-m) ...< n+1 : k+j*(k-m)<n-k+1 =>
% j<(n-2*k+1)/(k-m) j = floor((n-k+1)/(k-m))

    thr = cell(1,numLoops);
    btm = cell(1,numLoops);
    % main loop to compute the theories in batches of length k.
    for j = 0:numLoops-1 % this is k-m+1 when computing the correlation. Maybe want to keep the same notation, and have k-m+1 here too?
        % check if first or last index
%         j
        if j==0 
            if j== numLoops-1 % if it is only a single loop
                tempThry = [ones(extraL,1); ts.Data;ones(extraR,1)]; % add circ
            else
                tempThry = [ones(extraL,1); ts.Data(1:k-extraL)];
            end
            if circular
                tempThry(1:extraL) = ts.Data(end-extraL+1:end);
            end
        else % might be two cases when we loop over!
            if j==numLoops-1 || (j==numLoops-2 &&  k+j*(k-m)-extraL > length(ts.Data))
                numPts = min(extraR, k+j*(k-m)-extraL-length(ts.Data));
                tempThry = [ts.Data(j*(k-m)-extraL+1:end); ones(numPts,1)];
                if circular
                    tempThry(end-numPts+1:end) = ts.Data(1:numPts);
                end
            else
                tempThry = ts.Data(j*(k-m)-extraL+1:k+j*(k-m)-extraL); % theory shifted by m/2
            end
        end
        
        % add rand integers to the place of bad theory
        badThry = tempThry>4;
        tempThry(badThry) = randi(4,1,sum(badThry));
        
        % now compute tempThry
        x = compute_theory_wrapper(tempThry, sets);

        %The main trick of getting dot products in O(n log n) time
        % compute fft of x
        X = fft(x);

        % product
        Z = X.*Y;

        % inverser fourier
        z = ifft(Z);
        
        curT = zeros(1,length(cutPointsL{j+1}));
        curB =  zeros(1,length(cutPointsL{j+1}));
        % ts starts from j, and computes convolution from j+m/2 to j+k-1+m/2
        % cutPointsL has exact locations we need already
        for t = 1:length(cutPointsL{j+1})
            curT(t) = mean(z(cutPointsL{j+1}(t):cutPointsR{j+1}(t)));
            curB(t) = sum(badThry((cutPointsL{j+1}(t):cutPointsR{j+1}(t))-m/2));
        end
        thr{j+1} = curT;
        btm{j+1} = curB;
    end
    
    thry = cell2mat(thr);
    bitmask=cell2mat(btm);
end

