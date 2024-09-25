function [bgAll, bG, synthStr, synthStr2, theoryStruct,setsGen] = load_synth_data(N, pxpsf, mL, stdL, ...
    snvr, stdE, isC, tL, minL,seqID,sequences )
    %   load_synth_data - this script generates synthetic test example for
    %   bargrouping paper

    %   Args:
    %       N - number random fragments
    %       pxpsf - psf in pixels
    %       mL - mean length of fragment
    %       stdM - std of signal
    %       stdL - std of fragment length
    %       snvr - synthetic-noise variance ratio
    %       stdE - std extension factor
    %       isC - is circular
    %       tL - total length
    %       minL - minimum length
    
    %   Returns:



    % as an option, can pass theory barcode as last element, then
    % i.e. the pxpsf and tL does not matter.

%%
    rng("default")


    if nargin < 1
        N = 200;% number of random fragments
        pxpsf = 2.72; % psf
        mL = 850; % mean length of fragment
        stdL = 100; %std of length of fragment
%         snr = 0.5:0.5:5;
%         snvr = 0.2:0.2:2; % 
%         snr = 0.2;
        snvr = [0.18 0.2 0.25 0.3 0.4 0.5 0.75 1 2 4];
        stdE = 0.05; %0.02; % fragment length-rescale factor std 0.02
        isC = 1; % whether circular barcode
        tL = 10000; % rand length in px
        minL = 150; % minimum length / same as minoverlap?
%         seqID = 0;
%         sequences = [];
        sequences = dir('files/*.fasta');
        seqID = 1;
    end

    stdM = 20; % additive noise mean

    bgAll = [];
    %     
    bG = cell(1,length(snvr));
    synthStr = cell(1,length(snvr));
    synthStr2 = cell(1,length(snvr));
    theoryStruct = cell(1,length(snvr));

    import Core.gen_synth_data;
    import Core.synth_to_struct;


    for i=1:length(snvr)
    
        % signal to noise ratio defined as a ratio of variances of signal
        % and bg (bg has zero mean intensity)
        stdN = sqrt(stdM.^2/snvr(i)); % additive noise std. Signal mean=200 is hardcoded

        % synthetic (not based on theory)
        [bG{i}, barLong, origPos, origFlip, origStr] = gen_synth_data(tL, pxpsf, ...
            N, mL, stdL, stdM, stdN, stdE, isC,seqID,sequences );

        % generate an output structure similar to that for real experiments
        [synthStr{i},synthStr2{i},theoryStruct{i}] = synth_to_struct(bG{i}, barLong, origPos, origFlip, origStr,tL,minL);
%         synth_to_struct( origPos, origFlip,origStr) 
    end

    setsGen = struct('N',N,'pxpsf',pxpsf,'mL',mL,'stdM',stdM,'stdL',stdL,'snvr',snvr,'stdE',stdE,'isC',isC,'tL',tL,'minL',minL,'seqID',seqID,'sequences',sequences);

end

