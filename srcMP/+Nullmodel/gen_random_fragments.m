function [barcodeGen, barLong, origPos, origFlip, origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, ...
    NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD, IS_CIRCULAR, bar_temp)
%   import Microscopy.Simulate.Core.apply_point_spread_function

  % Args:
    %     TOTAL_RAND_LENGTH - total barcode length
    %     PSF_WIDTH_PIXELS - point spread function in pixels
    %     NUM_RAND_FRAGMENTS - number of random fragments
    %       MEAN_FRAGMENT_LENGTH - mean fragment length
    %       FRAGMENT_LENGTH_STD - fragment length std
    %     ADDED_NOISE_MEAN - noise mean
    %     ADDED_NOISE_STD - noise std  
    %     FRAGMENT_STRETCH_STD - fragment stretch std
    %     IS_CIRCULAR - whether barcode considered is circular
    %     BAR_LONG - input barcode, if already known.
  
    % Example
    % TOTAL_RAND_LENGTH = 10000;
    % PSF_WIDTH_PIXELS = 2.7;
    % NUM_RAND_FRAGMENTS = 100;
    % MEAN_FRAGMENT_LENGTH = 850;
    % FRAGMENT_LENGTH_STD = 100;
    % ADDED_NOISE_MEAN = .4;
    % ADDED_NOISE_STD = .1;
    % FRAGMENT_STRETCH_STD = 0.02;
    % sF = [0.95:0.01:1.05];
    % IS_CIRCULAR=1;


  % Return:
  %
  % fragments - fragments from the long barcode
  % barLong - long barcode
  % origPos - position on the long barcode, 
  % origFlip - barcode orientation o in {1,-1} - 1 - normal orientation, -1 -
  % reverse
  
  
    % initialize outputs
    barcodeGen = cell(NUM_RAND_FRAGMENTS, 1); 


    origPos = zeros(1, NUM_RAND_FRAGMENTS);
    origFlip = zeros(1, NUM_RAND_FRAGMENTS);
    origStr = zeros(1, NUM_RAND_FRAGMENTS);

    % create long barcode
    EXTRA_PX = round(6*PSF_WIDTH_PIXELS);
    
    if nargin < 10
        bar_temp = randn(1, TOTAL_RAND_LENGTH)+200;
    end
    
    if IS_CIRCULAR
%         barTemp = randn(1, TOTAL_RAND_LENGTH)+10;
        barLong  = [bar_temp bar_temp];
    else
        barLong = bar_temp;
    end
    MIN_LEN = 3*EXTRA_PX;

    
    % loop through number of fragments
    for i = 1:NUM_RAND_FRAGMENTS
        % length of short barcode. add 12*PSF since we'll want it linear
        barShortLen = max(MIN_LEN, round(MEAN_FRAGMENT_LENGTH + randn * FRAGMENT_LENGTH_STD));
    
        % position along the long barcode.  
        if IS_CIRCULAR
            barShortPos = randi(TOTAL_RAND_LENGTH);
        else
            barShortPos = randi(TOTAL_RAND_LENGTH - barShortLen + 1);
        end
        
        % now extract from the long barcode (together with 2*EXTRA_PX
        barShort = [zeros(1,EXTRA_PX) barLong(barShortPos:barShortPos + barShortLen - 1) zeros(1,EXTRA_PX)];
        % 
  
        % re-scale factor
        barShortStretch = 1 + randn * FRAGMENT_STRETCH_STD;
        % re-scaled length
        barShortStretchedLen = round(barShortStretch * barShortLen);
        
        % re-scale barcode
        barShort = interp1(barShort, linspace(1, barShortLen, barShortStretchedLen));
        
        % add noise
        barShort = barShort + randn(1, barShortStretchedLen) * (ADDED_NOISE_MEAN + randn * ADDED_NOISE_STD);
        
        % apply PSF
        barShort = imgaussfilt(barShort,PSF_WIDTH_PIXELS,'Padding','circular');

        barShort = barShort(EXTRA_PX+1:end-EXTRA_PX);
    
        origPos(i) = barShortPos;
        origStr(i) = barShortStretch;
        origFlip(i) = rand >= .5;

        if origFlip(i) 
          barcodeGen{i}.rawBarcode = fliplr(barShort);
        else
          barcodeGen{i}.rawBarcode = barShort;
        end
        barcodeGen{i}.rawBitmask = logical(zeros(1,length(barShort)));
        barcodeGen{i}.rawBitmask(round(EXTRA_PX/2)+1:end-round(EXTRA_PX/2))=1;
    end

    if IS_CIRCULAR
        barLong =  imgaussfilt(barLong(1:TOTAL_RAND_LENGTH),PSF_WIDTH_PIXELS,'Padding','circular');
    else
        barLong =  imgaussfilt([zeros(1,EXTRA_PX) barLong zeros(1,EXTRA_PX)],PSF_WIDTH_PIXELS,'Padding','circular');
        barLong = barLong(EXTRA_PX+1:end-EXTRA_PX);
    end
    

end
%   apply_point_spread_function(barLong, PSF_WIDTH_PIXELS, 1);
%   barLongr = barLong(round(6*PSF_WIDTH_PIXELS+1):round(end-6*PSF_WIDTH_PIXELS));
