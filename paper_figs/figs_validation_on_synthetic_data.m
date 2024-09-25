% validate ground-truth

% Does not really have much to do with the assembly, but validating the
% validation.

% Generate simulated barcodes for given SNR
% Compare to theory
% Calculate MP length based statistics
% Report how many end up at correct place

% 1) Validate the bargrouping_minimum_length procedure, which picks the
% thresh-length based on global sequence similarity. One can take a similar
% sequence (same species?), to estimate this. When running global
% alignment, this could be done for each position (since we find where it
% matches)


% Over different signal-to-noise-ratio

pathMain = '/home/avesta/albertas/reps/'; % path to reps.


% addpath(genpath([pathMain,'lldev']))
addpath(genpath([pathMain,'hca']))
addpath(genpath([pathMain,'bargroupingprototype']))
addpath(genpath([pathMain,'discriminative_sequences']))


timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% Generate synthetic data
N = 200;% number of random fragments
pxpsf = 2.72; % psf
mL = 850; % mean length of fragment
stdL = 100; %std of length of fragment
snr = 0.5:0.5:5;
stdE = 0.05;%0.02; % fragment length-rescale factor std 0.02
isC = 1; % whether circular barcode
tL = 10000; % rand length in px
minL = 150; % minimum length / same as minoverlap?

% synthetic data (if not generated yet
import Core.load_synth_data;
[bgAll, bG, synthStr, synthStr2, theoryStruct] = load_synth_data(N,2.72,850,100,snr,0.05,1,10000,150,0,[]);


%% Validation 1 - make sure that minLen criteria gives no/few false positives. FigS10.eps
import Validation.bml_validation;
[tpr, fpr] = bml_validation(synthStr,bG,theoryStruct,N, pxpsf, mL, stdL, ...
    snr, stdE, isC, tL, minL );



%% Validation 2 - validation overlap 
import Validation.ol_validation;
ol_validation(synthStr,theoryStruct, bG,sF,minOverlap)


%% Validation 3 - assembly
import Validation.a_validation;
a_validation() % might need some things from previous

% Validation 4 - reference based validation
import Validation.rb_validation;

rb_validation() % here we generate contigs based on consensus from reference (so we map to reference, etc.)

% validation 5 - find best overlap length
import Validation.length_Validation;
length_Validation