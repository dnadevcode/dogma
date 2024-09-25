function [barcodeGen, barLong, origPos, origFlip, origStr] = gen_synth_data(tL, pxpsf, ...
    N, mL, stdL, stdM, stdN, stdS, isC, seqID, sequence, bar_temp)

  % Args:
    %       N - number random fragments
    %       pxpsf - psf in pixels
    %       mL - mean length of fragment
    %       stdM - std of signal
    %       stdL - std of fragment length
    %       snr - signal to noise ratio
    %       stdE - std extension factor
    %       isC - is circular
    %       tL - total length
    %       minL - minimum length
  
    % Example
    %     N = 100;% number of random fragments
    %     pxpsf = 2.72; % psf
    %     mL = 850; % mean length of fragment
    %     stdL = 100; %std of length of fragment
    %     stdM = 2; % additive noise mean
    %     snr = 1:10;
    %     stdE = 0.05;%0.02; % fragment length-rescale factor std 0.02
    %     isC = 1; % whether circular barcode
    %     tL = 10000; % rand length in px
    %     minL = 150; % minimum length / same as minoverlap?


  % Return:
  %
  % fragments - fragments from the long barcode
  % barLong - long barcode
  % origPos - position on the long barcode, 
  % origFlip - barcode orientation o in {1,-1} - 1 - normal orientation, -1 -
  % reverse
  
  
    % initialize outputs
    barcodeGen = cell(N, 1); 


    origPos = zeros(1, N);
    origFlip = zeros(1, N);
    origStr = zeros(1, N);

    % create long barcode
    EXTRA_PX = round(6*pxpsf);
    
    if nargin < 12
        % evently distributed sequence
        [sets] = Core.Default.read_default_sets('hcasets.txt');
        sets.theoryGen = sets.default;
        import CBT.Hca.Core.Theory.choose_cb_model; % needs HCA
        [sets.model ] = choose_cb_model(sets.theoryGen.method,sets.theoryGen.pattern);

        sets.theoryGen.psfSigmaWidth_nm = 300;   % psf width in nm
        sets.theoryGen.pixelWidth_nm = 110;
        sets.theoryGen.meanBpExt_nm = 0.3;       % mean base-pair extension
        sets.theoryGen.k = 200000;                         %length of the small fragment to do fft on
        sets.theoryGen.m = 26000;                          % number of base-pairs to leave at the edges
        sets.theoryGen.isLinearTF = 0;               % is theory circular
        sets.theoryGen.addpsf = 0;               % is theory circular
        sets.theoryGen.computeBitmask = 0;
        if seqID==0
            name.Data = randi(4,round(tL*sets.theoryGen.pixelWidth_nm/sets.theoryGen.meanBpExt_nm),1); % even distribution
            name.Filename = [];
        else
            fasta =  fastaread(fullfile(sequence(seqID).folder,sequence(seqID).name));
            name.Data = nt2int(fasta.Sequence)';
            name.Filename = []; 
        end
        import CBT.Hca.Core.Theory.compute_theory_barcode; % needs HCA
        [ theory, header, bitmask] = compute_theory_barcode( name,sets);
        bar_temp = stdM*zscore(theory)+200;

%         bar_temp = stdM*randn(1, tL)+200;
    end

    tL = length(theory);
    
    if isC
%         barTemp = randn(1, TOTAL_RAND_LENGTH)+10;
        barLong  = [bar_temp bar_temp];
    else
        barLong = bar_temp;
    end
    MIN_LEN = 3*EXTRA_PX;



    
    % loop through number of fragments
    for i = 1:N
        % length of short barcode. add 12*PSF since we'll want it linear
        barShortLen = max(MIN_LEN, round(mL + randn * stdL));
    
        % position along the long barcode.  
        if isC
            barShortPos = randi(tL);
        else
            barShortPos = randi(tL - barShortLen + 1);
        end
        
        % now extract from the long barcode (together with 2*EXTRA_PX
        barShort = [zeros(1,EXTRA_PX) barLong(barShortPos:barShortPos + barShortLen - 1) zeros(1,EXTRA_PX)];
        % 
  
        % re-scale factor
        barShortStretch = 1 + randn * stdS;
        % re-scaled length
        barShortStretchedLen = round(barShortStretch * barShortLen);
        
        % re-scale barcode
        barShort = interp1(barShort, linspace(1, barShortLen, barShortStretchedLen));
        
        % add noise
        barShort = barShort + stdN*randn(1, barShortStretchedLen);% * (0 + randn * stdN); % noise has zero mean (i.e. perhaps we have removed noise)
        
        % apply PSF
        barShort = imgaussfilt(barShort,pxpsf,'Padding','circular');


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

    if isC
        barLong =  imgaussfilt(barLong(1:tL),pxpsf,'Padding','circular');
    else
        barLong =  imgaussfilt([zeros(1,EXTRA_PX) barLong zeros(1,EXTRA_PX)],pxpsf,'Padding','circular');
        barLong = barLong(EXTRA_PX+1:end-EXTRA_PX);
    end
    
end


        %% Analysis of SNR
%         st = [];
%         newVar = [];
%         for i=1:10000
%         barN = stdN*randn(1, barShortStretchedLen);
%         barN2 = imgaussfilt(barN,pxpsf,'Padding','circular');
%         st = [st var(barN2(EXTRA_PX:end-EXTRA_PX))];
% 
% %         st = [st std(barN(EXTRA_PX:end-EXTRA_PX))/std(barN2(EXTRA_PX:end-EXTRA_PX))];
% %         newVar = [newVar var(barN)+var(barN2)];
%         end
% 
%         stlong = imgaussfilt(barLong,pxpsf,'Padding','circular');
%         var(stlong)/mean(st)
% 
% %         [mean(st) stdN.^2/sqrt(1+((stdN.^2)/(pxpsf.^2)))]
% 
%         figure,histogram(st)
%             %