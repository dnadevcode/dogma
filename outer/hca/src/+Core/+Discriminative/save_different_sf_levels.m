function [] = save_different_sf_levels(inputData,inputTable)

% % Input data for all coefficients and input PCC table for all coefficients
% inputData = '20221031_Sample 398-st1_1463.1bpPERpx_0.170nmPERbp0_0.9sf_allcoefs.mat';
% inputTable = '20221031_Sample 398-st1_1463.1bpPERpx_0.170nmPERbp_MP_w=0_table_2024-04-29_17_11_45.txt';

if nargin < 2
    % todo: these two files should be generated from the same run, so name
    % them the same in the future (add time to allcoefs file too)
    [inputData,inputFold] = uigetfile('*.mat','Select allcoefs file');
    [inputTable,inputFold2] = uigetfile('*.txt','Select PCC file',inputFold);

    inputTable = fullfile(inputFold2,inputTable);
    inputData = fullfile(inputFold,inputData);

end

import Core.export_coefs_local;

outdir = '';

load(inputData);

tic
ds = datastore(inputTable);
tt = tall(ds);
thryNames = gather(tt.Var1);
bnames = ds.SelectedVariableNames(2:4:end);
toc
inputstr = strsplit(inputData,'MP_');


% import Core.load_local_alignment_results_from_files;
% [~,bnames,~,thryNames] = load_local_alignment_results_from_files(foldername);

maxCoef = cell(1,size(matAllCoefs,1));
maxOr = cell(1,size(matAllCoefs,1));
maxPos = cell(1,size(matAllCoefs,1));
maxSecondPos = cell(1,size(matAllCoefs,1));
maxlen = cell(1,size(matAllCoefs,1));
bestSF = cell(1,size(matAllCoefs,1));
for i=2:size(matAllCoefs,2)/2 % 1 would be original level
    %% Now find the best coefficient from matAllCoefs (using a cascading or whatever scheme)
    sfLevel = i;
    for barid =1:size(matAllCoefs,1)
        [singleCoef , singlePos ] =  max(matAllCoefs(barid,sfLevel:end-sfLevel+1,:),[],2);
        pos  = squeeze(singlePos)';
        maxCoef{barid} =  squeeze(singleCoef);
        maxOr{barid} = zeros(size(maxCoef{barid}))';
        maxPos{barid} = zeros(size(maxCoef{barid}))';
        maxSecondPos{barid} = zeros(size(maxCoef{barid}))';
        maxlen{barid} = zeros(size(maxCoef{barid}))';
        bestSF{barid} =  zeros(size(maxCoef{barid}))';
    end
%     toc



export_coefs_local(thryNames',maxCoef,maxOr,maxPos,maxlen, bestSF, bnames,[outdir,inputstr{1}, '_MP_w=',num2str(sets.w),'_','sfDepth_',num2str(i),'_']);


end

end

