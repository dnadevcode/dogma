function [theory, Y, cutPointsL, cutPointsR, extraL, extraR, ...
    psfSigmaWidth_bps, pxSize, k, m, circular, n, nSeq, numPx,curIdx, leftOverTheory,bpPlaceL,bpPlaceR] =...
    set_up_theory_for_calculation_gc(sets, sequence)
    %
    % set_up_theory_for_calculation
    %
    %   Args:
    %       sets : settings with theory parameters
    %       sequence : dna sequence used (or could be memfile, if memory mapped)
    %
    %
    %   Returns:
    %
    %
    % convert nt  to int
    
    % Convert sequence to int, this can be already done if data is
    % mem-mapped
%     sequence = nt2int(sequence)';


    %% compute CG ratio - this will be used to calculate how many bp to take 
    % for each pixel.
%     numNT = numel(sequence.Data);
%     blockSize = 500; % how many bp to make the bloc size?
%     numBlocks = floor(numNT/blockSize);
% 
%     ratio = zeros(numBlocks+1,1);
% 
%     A = nt2int('A'); C = nt2int('C'); G = nt2int('G'); T = nt2int('T');
% 
%     for count = 1:numBlocks
%         % calculate the indices for the block
%         start = 1 + blockSize*(count-1);
%         stop = blockSize*count;
%         % extract the block
%         block = sequence.Data(start:stop);
%         % find the GC and AT content
%         gc = (sum(block == G | block == C));
%         at = (sum(block == A | block == T));
%         % calculate the ratio of GC to the total known nucleotides
%         ratio(count) = gc/(gc+at);
%     end
% 
%     block = sequence.Data(stop+1:end);
%     gc = (sum(block == G | block == C));
%     at = (sum(block == A | block == T));
%     ratio(end) = gc/(gc+at);
% 
%     strF = 0.1;
%     dFactors = (ratio/0.5-1)*strF+1;
%     
%     dFactors(isnan(dFactors)) = 1;
%     dFactors = dFactors/sum(dFactors)*(length(ratio));
    
    % now create a vec where 
    
    % each dFactor describes the re-scaling factor for blockSize
    % base-pairs.
    % Example
    % if dFactors=1.1, this means that blockSize base-pairs correspond to 
    % 1.1 blockSize base-pairs.
    % since each dFactors is representing blockSize, in the end we
    % normalize that the total number of blockSize would be the same
    
    % edge cases are when we take a pixel value between the two blocks

%     ts = sequence;
    %% Parameters used
    % psf in basepairs, we need to convolve with such Gaussian
    psfSigmaWidth_bps = sets.theoryGen.psfSigmaWidth_nm / sets.theoryGen.meanBpExt_nm;
    % This function converts bpRes to pxRes
    pxSize = sets.theoryGen.pixelWidth_nm/sets.theoryGen.meanBpExt_nm;
    % maybe no need to re-write these?    
    k = sets.theoryGen.k; % length of the small fragment to do fft on
    m = sets.theoryGen.m; % number of base-pairs to leave at the edges
    circular = ~sets.theoryGen.isLinearTF;
    
    nSeq = length(sequence.Data); % theory + extra things on the left and right

%     nSeq+m+2*round(pxSize/2)+l = k+(j+1)*(k-m), l<k-m
%  old l =  k+j*(k-m)-nSeq+m/2+round(pxSize)+m/2+round(pxSize)
%  nSeq+m+l = k+j*(k-m) => l =  k+j*(k-m)-nSeq-m
    % j:  k+j*(k-m) <= nSeq=> j <= (nSeq-k)/(k-m)
    j =  ceil((nSeq-k+m)/(k-m)); % number of j steps (starts at 0)
    l = k+(j)*(k-m)-nSeq-m;%-m; % add extra to have a single loop
    
    extraL = m/2; %+ceil(pxSize); % extra bp to add at idx=1 
    extraR = m/2+l; %+ceil(pxSize/2)+l; % extra bp to add at idx=j
    
	n = extraL+nSeq+extraR; % total number of bp considered

    numPx = ceil(nSeq/pxSize); % number of pixel
    
    % initialize theory barcode
    theory = zeros(1,numPx);
   
    cutPointsL = cell(1,j+1);
    cutPointsR = cell(1,j+1);
    
    blockSize=sets.blockSize;
    strF=sets.strF;
    %%
    % pixel positions (on extended theory:
    % m/2+1+ceil(pxSize/2) (first pixel), m/2+ceil(pxSize)+nSeq+ceil(pxSize/2) (last
    % pixel), and we don't allow the last px to overflow 
    bpPlaceL = [];
    bpPlaceR = [];
    for jIdx=0:j
        % in this frame the pixels are labelled 1:k, but they represent 
        %   1:k=>      j*(k-m)+1:k+j*(k-m)
        % interesting pixels: m/2+1:k-m/2 (so k-m)

        % these start and stop are as if on the real theory
        start = jIdx*(k-m)+1; % ignore the m/2 in the beginning for the moment
        if jIdx==j
            % in this case
            stop = nSeq;%-ceil(pxSize/2); % we maybe discard the last unfilled pixel
        else
            stop = k+jIdx*(k-m)-m; % stop, so k-m bp are taken
        end
        
%         [stop start]
        
%         if jIdx==0
% %             dFactors = calc_gc(sequence.Data(start-extraL:stop),500,0.1);
%         else
% blockSize=5000; % dependence on blocksize
        [dFactors] = calc_gc(sequence.Data(start:stop),blockSize,strF);
%         end
%         dFactors = ones(66,1)
        
        % TODO: instead check based on previous stop position
        % correct strange behaviour that the first bp actually starts at
        % -ceil(pxSize/2):
        startPos = floor((start+ceil(pxSize/2))/pxSize); %-1; % start pixel
        stopPos = floor((stop-ceil(pxSize/2))/pxSize); % stop pixel
        
             
        % first loop is to calculate newpx factors, i.e. new pxSizes
        pxSizes = zeros(1,stopPos-startPos+1);
        pxSizes(1) = pxSize*dFactors(1);
        pxPosNext = -pxSize/2+1+pxSizes(1);
        for kIdx=2:stopPos-startPos+1
            pxSizes(kIdx) = pxSize*dFactors(max(1,(floor(pxPosNext/blockSize)+1))); % pixel has less basepairs if gc content is higher
            pxPosNext = pxPosNext+pxSizes(kIdx);
        end
        
        dFactors = dFactors*pxSize/mean(pxSizes);
        
%          % first loop is to calculate newpx factors, i.e. new pxSizes
%         pxSizes = zeros(1,stopPos-startPos+1);
%         pxSizes(1) = pxSize*dFactors(1);
%         pxPosNext = -pxSize/2+1+pxSizes(1);
%         for kIdx=2:stopPos-startPos+1
%             pxSizes(kIdx) = pxSize*dFactors(max(1,(floor(pxPosNext/blockSize)+1))); % pixel has less basepairs if gc content is higher
%             pxPosNext = pxPosNext+pxSizes(kIdx);
%         end
%         
%         dFactors = dFactors*pxSize/mean(pxSizes);
% 
%           pxSizes = zeros(1,stopPos-startPos+1);
%         pxSizes(1) = pxSize*dFactors(1);
%         pxPosNext = -pxSize/2+1+pxSizes(1);
%         for kIdx=2:stopPos-startPos+1
%             pxSizes(kIdx) = pxSize*dFactors(max(1,(floor(pxPosNext/blockSize)+1))); % pixel has less basepairs if gc content is higher
%             pxPosNext = pxPosNext+pxSizes(kIdx);
%         end
%         
%         dFactors = dFactors*pxSize/mean(pxSizes);

       
%            pxSizes = zeros(1,stopPos-startPos+1);
%         pxSizes(1) = pxSize*dFactors(1);
%         pxPosNext = -pxSize/2+1+pxSizes(1);
%         for kIdx=2:stopPos-startPos+1
%             pxSizes(kIdx) = pxSize*dFactors(max(1,(floor(pxPosNext/blockSize)+1))); % pixel has less basepairs if gc content is higher
%             pxPosNext = pxPosNext+pxSizes(kIdx);
%         end
%         
%         dFactors = dFactors*pxSize/mean(pxSizes);
%         
%         
% 
%        
%            pxSizes = zeros(1,stopPos-startPos+1);
%         pxSizes(1) = pxSize*dFactors(1);
%         pxPosNext = -pxSize/2+1+pxSizes(1);
%         for kIdx=2:stopPos-startPos+1
%             pxSizes(kIdx) = pxSize*dFactors(max(1,(floor(pxPosNext/blockSize)+1))); % pixel has less basepairs if gc content is higher
%             pxPosNext = pxPosNext+pxSizes(kIdx);
%         end
%         
%         dFactors = dFactors*pxSize/mean(pxSizes);
%         % to check the pxSizes
%         % second loop with corrected dFactors to have 
%         pxSizes = zeros(1,stopPos-startPos+1);
%         pxSizes(1) = pxSize*dFactors(1);
%         pxPos = zeros(1,stopPos-startPos+1);
%         pxPos(1) = startPos*pxSize-pxSize/2 +1;
%         pxPosNext = -pxSize/2+pxSizes(1);
%         for kIdx=2:stopPos-startPos+1
%             pxSizes(kIdx) = pxSize*dFactors(max(1,(floor(pxPosNext/blockSize)+1))); % pixel has less basepairs if gc content is higher
%             pxPosNext = pxPosNext+pxSizes(kIdx);
%         end
% %         
        
        
        % curStrFactor
%         -(jIdx*(k-m))+m
      
        % compute gc content for thes
        pxSizeCur = pxSize*dFactors(1); % pixel has less basepairs if gc content is higher
        bpPlaceL{jIdx+1}(1) = start-pxSize/2+1;
        bpPlaceR{jIdx+1}(1) = start-pxSize/2+pxSizeCur;
        px = zeros(1,stopPos-startPos+1);
        px(1) = pxSizeCur;
        % redo bpPlaceL and bpPlaceR:
        for kIdx=2:stopPos-startPos
%             pxSizeCurrent = pxSize*(dFactors(kIdx
            bpPlaceL{jIdx+1}(kIdx) = bpPlaceR{jIdx+1}(kIdx-1)+1 ; %or pxSize
            pxSizeCur = pxSize*dFactors(max(1,(floor((bpPlaceL{jIdx+1}(kIdx)-start)/blockSize)+1))); % pixel has less basepairs if gc content is higher
%             pxSizeCur
            px(kIdx) = pxSizeCur;
            bpPlaceR{jIdx+1}(kIdx) =  bpPlaceL{jIdx+1}(kIdx)+pxSizeCur-1;   
        end
        bpPlaceL{jIdx+1}(end) =  bpPlaceR{jIdx+1}(end-1)+1;
        bpPlaceR{jIdx+1}(end) = stop-pxSize/2;
%         bpPlaceL{jIdx+1} =  floor((startPos:stopPos)*pxSize-pxSize/2+1); %or pxSize
%         bpPlaceR{jIdx+1} =  floor((startPos:stopPos)*pxSize+pxSize/2);   
%         mean(pxSizeCur)

        
%         % PX:
%         startPos = floor((start+ceil(pxSize/2))/pxSize); %-1; % start pixel
%         stopPos = floor((stop-ceil(pxSize/2))/pxSize); % stop pixel
%         
%         bpPlaceL{jIdx+1} =  floor((startPos:stopPos)*pxSize-pxSize+1); %or pxSize
%         bpPlaceR{jIdx+1} =  floor((startPos:stopPos)*pxSize+pxSize);   
%         
        
      
        % start position on theory extract:
        cutPointsL{jIdx+1} =  floor(bpPlaceL{jIdx+1}-(jIdx*(k-m))+m);%+m gives a position
        cutPointsR{jIdx+1} =  min(k,floor(bpPlaceR{jIdx+1}-(jIdx*(k-m))+m));   
        % todo: check again with the convolutions
        
    end
%     figure,plot(bpPlaceR{1}-bpPlaceL{1})
%     
%     floor(bpPlaceL{2}(1))
%     floor(bpPlaceR{1}(end))
%     floor(bpPlaceR{2}(end))
    
%     floor(bpPlaceL{2003}(1))
%     floor(bpPlaceR{2002}(end))

%     
%     floor(bpPlaceL{1}(1))
%     floor(bpPlaceR{1}(1))
%         
%     floor(bpPlaceL{1}(end))
%     floor(bpPlaceR{1}(end))
%     length(bpPlaceL{1})
    %%
%         % start: m/2+ceil(pxSize/2) is the firs
% %         start = jIdx*(k-m)+1; % position of center pixel  extraL..k (so k-m/2-ceil(pxSize/2) last possible bp we take for pixel which starts at extraL+1)
%         % so we check pixels in extraL... k-m/2-ceil(pxSize/2) 
%         %         jIdx*(k-m)+1+extraL
% % % %          if jIdx==j
% % % %             stop = nSeq+extraL+ceil(pxSize/2); % stop  k+(j+1)*(k-m)
% % % %          else
% %             stop = k+jIdx*(k-m)-extraL-m/2-ceil(pxSize/2); % stop
% %          end
%         jIdx=0;
%         jIdx*(k-m)+1
%         k+jIdx*(k-m)-m/2-ceil(pxSize/2)
%         jIdx=1;
%         jIdx*(k-m)+1
%         k+jIdx*(k-m)-extraL-m/2-ceil(pxSize/2)
%          % globally? 
%         startPos = floor(start/pxSize); %-1; % start pixel
%         stopPos = floor(stop/pxSize); % stop pixel
% 
%         bpPlaceL{jIdx+1} =  floor((startPos:stopPos)*pxSize-pxSize/2+1); %or pxSize
%         bpPlaceR{jIdx+1} =  floor((startPos:stopPos)*pxSize+pxSize/2);   
%         

%     end
        % and which bp values from the map z(..) we have to take
%         
%         % we take j*(k-m)+1:k+j*(k-m). These represent
%     %     j*(k-m)+1+m/2: k+j*(k-m) - m/2
%         stop = (jIdx+1)*(k-m);
%         stopPx = floor(stop/pxSize);
%     %     % test for start / should only check edge cases to see if same pixel
%     %     can be mapped twice sometimes
%     %     floor( ceil(start/pxSize)*pxSize-pxSize+1)
%     %     % test for stop
%     %     floor(stop/pxSize)*pxSize
%  
%     end
    % last block will not be full!
    
    % note that  cutPointsL(i):cutPointsR(i) will not always be the same
    % and changes a little bit, because we take the floor value ( so can
    % differ in +-1 basepair)
    
%     cutPoints = [1 2*floor(pxSize):floor(pxSize):nSeq]; % how best to deal with this. Don't want to have different
        
    % current index for px theory
    curIdx = 1;
        
    % initialize kernel (of length k)
    y = zeros(k,1);
    % kernel values
    y(1:m) = images.internal.createGaussianKernel(psfSigmaWidth_bps, m);

    %    Reverse query: not needed if kernel is symetric!
    %     y = ker(end:-1:1); %Reverse the query
    %     y(m+1:k) = 0; %append zeros

    % Fourier transform of the kernel. 
    Y = fft(y);
        
    % Only save the pixelated theory values in vector "theory"

    leftOverTheory = [];
    
    
    
end


    % total length:
%     jMax = (j+1)*(k-m);% floor((n-2*k+1)/(k-m))+1; % number of loops we divide theory calculation into

    % on the left we can add as much as we want. We add m/2 since the
    % Gaussian we use will have half length m/2. Add ceil(pxSize) since the
    % first pixel starts at m/2+ceil(pxSize/2), but later pixels can start
    % anywhere between m/2 and m/2+ceil(pxSize) - so we make sure that full
    % Gaussian is always used
    
    % Sometimes sequence can be very long, we don't want to copy it to a
    % new sequence. Instead we'll just take specific values on the left and
    % the right..
 %     ts = [ones(extraL,1); sequence; ones(extraR,1)];
%     if circular % so if circular, we add the theory data instead of simple ones
%         ts(1:extraL) = sequence(end-extraL+1:end);
%         ts(end-extraR+1:end) =  sequence(1:extraR);
%     end
        
    % This is moved to outer structure avoid repetitive (periodic) things in the theory barcode.
    %     s = rng;
    %     rng(0,'twister');
    %     ts(ts > 4) = randi(4,1,sum(ts > 4));
    %     rng(s);
    
    
    
    % follow notation from convert_bpRes_to_pxRes
%     theory(floor(nSeq/pxSize)) = 0;
    % cutPoints. after reaching these points, the value in a pixel of
    % theory is computed. the first one starts from 1. We don't end exactly
    % and nSeq, so there's some information loss (almost up to a pixel) 
    % especially if the barcode is circular. 
 
    % find the cutpoints in theory at which should we compute the theory
%     cutPointsR = ones(1,numPx);
%     cutPointsL = ones(1,numPx);
%     for i = 1:length(cutPointsR)
%         cutPointsR(i) = floor((i*pxSize)+pxSize); % so we take [-pxSize:pxSize] around the value of interest
%         cutPointsL(i) = floor((i*pxSize)-pxSize+1);
%     end
    
%     jMax =  floor((n-2*k+1)/(k-m))+1; % number of loops we divide theory calculation into


%     left = floor((0:1000)*pxSize-pxSize/2+1);
%     right= floor((0:1000)*pxSize+pxSize/2);



        % these are mapped to m:k-1 on the convolution vector. In practice
        % we want to have all pixels mapped, so we drop the last one and
        % add one to the left
        % 
        % when convolving with Gaussian, something that is placed at 1 has
        % a peak centered around m/2 in the convolution, so m/2+1 has a peak at m+1 - that's
        % the first bp of interest, and the last bp of interest is k-m/2-1,
        % which has a peak at k-1, therefore z(m:k-1)
        %
        % in jIdx=1..j+1 part it's the same. Just relation differs. We
        % allow the start and end to differ a little bit so that pixel
        % values would be taken from independent runs
        % 
        
        

                % last one should be just outside nSeq
        %         % importantly, which bp values for these
        %         floor(stPx*pxSize-pxSize/2+1)
        %         floor(stPx*pxSize+pxSize/2)
        %         floor(stopPx*pxSize-pxSize/2+1)
        %         floor(stopPx*pxSize+pxSize/2)
        %         if jIdx==j+1
        %             % in this case only 
        %             
        %         else
            %         end
      