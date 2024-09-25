function [bitmask, posY,mat,threshval,threshstd,badMol,bitWithGaps] = median_filt(km, filterM,bgSigma,threshval,threshstd, N,Nzero)

% use median filter to detect edges
%   Args:
%       km - cell of input kymograph with the same intensity settings
%       filterM - filter size for filtering
%
%   Returns:
%       posY - positions left/right
%       mat - output thresholded image
%       bitmask - bitmask of movie
%       threshval - mean bg
%       threshstd - std bg
%       badMol - which mols has bad mask

    if nargin < 2
        filterM = [5 15];
    end

    if nargin < 6 || isempty(N)
        N = 50;
    end

    if nargin < 7
        Nzero = 0;
    end
    
    if nargin < 4

        vecvals = cell(1,length(km));
        for i=1:length(km)
            vecvals{i} = km{i}(:);
        end
        fullmat = double(cell2mat(vecvals'));
        fullmat(fullmat< 1) = nan; % assume nan's everything that is 0
        fullmat = fullmat(~isnan(fullmat)); % remove nan's

        vals = kmeans(fullmat,2);
        % alternative: thresh based on min point between first two peaks

        % check which is larger
        mat1 = fullmat(vals==1);
        mat2 = fullmat(vals==2);

        if mean(mat1)<mean(mat2)
            bg = mat1;
        else
            bg = mat2;
        end

        % filter, same number of points as signal filter
        filterS = [1 filterM(1)*filterM(2)];

        % calculate mean and standard deviation
        threshval = nanmean(medfilt2(bg,filterS,'symmetric'));
        threshstd = nanstd(medfilt2(bg,filterS,'symmetric'));
    end
    % number of standard deviations to allow
    if nargin < 3
        bgSigma = 1;
    end
    
    badMol = zeros(1,length(km));
    mat = cell(1,length(km));
    bitmask = cell(1,length(km));
    bitWithGaps = cell(1,length(km));
    posY = cell(1,length(km));

    for i=1:length(km)
        
        
        K = medfilt2(km{i},filterM,'symmetric') > threshval+bgSigma*threshstd;
        
        % potential speed up : run medfilt2 for discretized image, x5 times
        % speed up
%         tic
%         discThresh = gray2ind((threshval+bgSigma*threshstd - min(km{i}))/(max(km{i})-min(km{i})),2^8);
%         tempKm = gray2ind((km{i}-min(km{i}))./(max(km{i})-min(km{i})),2^8);
%         K = medfilt2(tempKm,filterM,'symmetric') > discThresh;
%         toc

      

    %     K = medfilt2(km{i},filterM,'symmetric') > meanD+bgSigma*varD; % should give correct also for unmedian filtered. Sigma bigger though?

%         figure;    imshowpair(imresize(K,[200 500]),imresize(km{i},[200 500]), 'ColorChannels','red-cyan'  )
    %     f

    %     figure,imagesc(km{i})
        mat{i} = K;

        if nargout >= 3
        [labeledImage, numBlobs] = bwlabel(K);
        numZeroRows = sum(0==sum(labeledImage,2));
%         numZeroCols = sum(0==sum(labeledImage,1));

        props = regionprops(labeledImage, 'Area');
        sortedVals = sort([props.Area],'desc');


        if numZeroRows==0 &&(numBlobs == 1 || (length(sortedVals)>1 && sortedVals(2) < N))

            % if there are two regions, could split into two molecules.
            [maxArea, largestIndex] = max([props.Area]);
            labK = labeledImage==largestIndex; % either just max or create a loop here
            
            if sum(0==sum(labK,2))~= 0;
                badMol(i)=1;
                continue;
            end
            

            posY{i}.leftEdgeIdxs = arrayfun(@(x) find(labK(x,:) >0,1,'first'),1:size(labK,1));
            posY{i}.rightEdgeIdxs = arrayfun(@(x) find(labK(x,:) >0,1,'last'),1:size(labK,1)); 
            

            numZeroPx = sum(arrayfun(@(x) sum(labK(x, posY{i}.leftEdgeIdxs(x): posY{i}.rightEdgeIdxs(x))==0) ,  1:length(posY{i}.leftEdgeIdxs)));
            
            if numZeroPx > Nzero % remove barcodes which are fragmented
                posY{i} = [];
                badMol(i)=1;
            else
                
                bitmask{i} = zeros(size(K));
                bitWithGaps{i} = labK;

                for j=1:size(bitmask{i},1)
                    bitmask{i}(j,posY{i}.leftEdgeIdxs(j):posY{i}.rightEdgeIdxs(j)) = 1;
                end
            end
        else
            posY{i} = [];
            badMol(i)=1;
        end
        end

    end
    
% ix = 2;
% kymo = km{ix};




end

