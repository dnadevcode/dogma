%Generates and saves result figures, also generates some of supplementary
%figures. Requires .mat files with pregenerated data

load('synth_2024-09-05_09_03_37resRun.mat'); % nc_000913.3


% requires
    % Signal Processing Toolbox
    % Image Processing Toolbox
    % Statistics and Machine Learning Toolbox
    % Computer Vision Toolbox
    %  Deep Learning Toolbox
    %  Symbolic Math Toolbox
    %  Bioinformatics Toolbox

%% Fig3: Plot with graphs of barcode islands
sets.minOverlap = 300; % minimum overlap
snrv = setsGen.snvr;

final_figure_3(sets, resRun, bG,snrv)

%% Fig4
final_figure_4(sets, resRun, bG,synthStr,theoryStruct)

%% Fig5
final_figure_5(resRun,bG,theoryStruct,synthStr)

%% Fig6
% load('synth_2024-09-05_09_03_37resRun.mat'); % nc_000913.3
load('experiment_data_oS.mat');
load('experiment_data_barcodeGen.mat');

final_figure_6

%% Fig7
final_figure_7
