function [t,outputLocs] = run_pipeline_scores(dirName, selIdxs, depth, windowWidths, sF, timeFramesNr, thryFiles, savedir, thryIdx,preloadedTheory,preloaded_barGen)

% run_pipeline_scores - optimized function to run experiment vs theory
% comparisons

%   Args:
%       selIdxs(1) - folder index
%       selIdxs(2) - subfolder index
%       dirName
%       windowWidths - alignment window widths, 0 is full
%       sF - length re-scaling factors
%       timeFramesNr - number of timeframes, 0 is all
%       thryFiles - theory files

t0 = tic; % start timer

% assign to ix and iy
ix = selIdxs(1);
iy = selIdxs(2);

% pre-assign things if they were not defined
if nargin < 5 || isempty(windowWidths)
    windowWidths = 250:100:600;
end

if nargin < 4
    depth = 1;
end

if nargin < 2
    dirName = '/proj/snic2022-5-384/users/x_albdv/data/bargrouping/local/';
end

addpath(genpath(dirName)); % add to path the location of kymos

[~,~] =mkdir('output'); % create output directory


if nargin < 1
    ix = 1;
end

dirStruct = dir(dirName);
dirStruct(~[dirStruct.isdir]) = [];  %remove non-directories
dirStruct(ismember( {dirStruct.name}, {'.', '..'})) = [];  %remove . and ..


if depth == 1
    subDir = dir(fullfile(dirStruct(ix).folder,dirStruct(ix).name));
    subDir(ismember( {subDir.name}, {'.', '..'})) = [];  %remove . and ..
    subDir(find(~cell2mat({subDir.isdir}))) = [];
else
    subDir = dirStruct(ix);
end

%
% iy = 1; % most likely single run
sets.dirName = fullfile(subDir(iy).folder,subDir(iy).name);

if nargin < 8
    savedir = sets.dirName;
end

% smart extract nmPerbp from folder name
[sets.nmbp, nmpx] = Input.extract_extension_params(sets.dirName);

% spltName = strsplit(sets.dirName ,'_');
% spltName2 = strsplit(spltName{end},'nm');
% spltName3 = strsplit(spltName{end-1},'nm');
% sets.nmbp = str2double(spltName2{1});
% nmpx = str2double(spltName3{1});


% load theory
thryFileIdx = find(arrayfun(@(x) ~isempty(strfind(thryFiles(x).name,num2str(nmpx))),1:length(thryFiles)));
sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);

%%

files = dir(fullfile(sets.dirName,'kymos','*.tif'));

sets.kymosets.filenames = arrayfun(@(x) files(x).name,1:length(files),'un',false);
sets.kymosets.kymofilefold = arrayfun(@(x) files(x).folder,1:length(files),'un',false);

% simdata = 0;
sets.output.matDirpath = 'output';
sets.filterSettings.filter = 0;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 3*300/nmpx;

sets.minLen = 150;
sets.genConsensus  = 0;

%  following "Strain-level bacterial typing directly from patient
% samples using optical DNA mapping"

if nargin < 6
    sets.timeFramesNr = 20;
else
    sets.timeFramesNr = timeFramesNr;
end

if nargin < 5
    sets.theory.stretchFactors = 0.8:0.025:1; %as per 
else
    sets.theory.stretchFactors = sF;
end

sets.alignMethod = 1; % nralign. Can fail on some images?
sets.edgeDetectionSettings.method = 'Otsu';

%% START CALCULATION
% load kymograph data
if nargin >=11
    barGen = preloaded_barGen;
else
    import Core.load_kymo_data;
    [kymoStructs,barGen] = load_kymo_data(sets);
    
    % rescale
    import Core.rescale_barcode_data;
    [barGen] = rescale_barcode_data(barGen,sets.theory.stretchFactors);
end

% save(['bars.mat'],'barGen','kymoStructs','sets');

% figure,tiledlayout(ceil(sqrt(length(kymoStructs))),ceil(length(kymoStructs)/sqrt(length(kymoStructs))),'TileSpacing','none','Padding','none')
% for i=1:length(kymoStructs)
%     nexttile;        imshowpair(imresize(kymoStructs{i}.alignedMask,[200 500]),imresize(kymoStructs{i}.alignedKymo,[200 500]), 'ColorChannels','red-cyan'  );    title(num2str(i));
% % imshowpair(imresize(kymoStructs{i}.unalignedBitmask,[200 500]),imresize(kymoStructs{i}.unalignedKymo,[200 500]), 'ColorChannels','red-cyan'  )
% end

clear kymoStructs; % releases memory used up by kymoStructs
% tic

% Load theory. Takes
% load(thryFile); 
sets.theoryFile{1} = sets.thryFile;
sets.theoryFileFold{1} = '';
sets.theory.precision = 5;
sets.theory.theoryDontSaveTxts = 1;
if nargin >=10 && ~isempty(preloadedTheory)
    theoryStruct = preloadedTheory;
else
    import CBT.Hca.UI.Helper.load_theory;
    theoryStruct = load_theory(sets);
end

% extract from name
sets.theory.nmbp = sets.nmbp*mean(sets.theory.stretchFactors);

if nargin>=9
    theoryStruct = theoryStruct(thryIdx);
end

% 
% tic
% ticBytes(gcp);

import CBT.Hca.Core.Analysis.convert_nm_ratio;
theoryStruct = convert_nm_ratio(sets.theory.nmbp, theoryStruct, sets );




thryNames = {theoryStruct.name};
barcodeNames = cellfun(@(x) x.name,barGen,'un',false);
% tocBytes(gcp);
% toc

% convert thry to correct nm/px

%

% sets.w = 300;
% sets.comparisonMethod = 'mass_pcc';
sets.genConsensus = 0;
sets.filterSettings.filter = 0;
% globalov = 0;

% import Core.extract_species_name; % find e-coli
% [speciesLevel, idc] = extract_species_name(theoryStruct);


% import CBT.Hca.Export.export_cc_vals_table;
% [T] = export_cc_vals_table( theoryStruct, comparisonStructAll, barcodeGenC,matDirpath);


%%

% theoryStruct([refNumsMP{5}]).name;
% windowWidths = 400:100:600;
sets.comparisonMethod = 'mpnan';


import CBT.Hca.Core.Comparison.hca_compare_distance;
import Core.export_coefs_local;

dateVec =  datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

outputLocs = cell(1,length(windowWidths));

for wIdx = 1:length(windowWidths)
    sets.w = windowWidths(wIdx);
    display(['Running w = ',num2str(sets.w)]);
    % only local lengths for which all length re-scaled versions passes the
    % threshold
    passingThreshBars = find(cellfun(@(x) sum(x.rawBitmask),barGen)*sets.theory.stretchFactors(1) >= sets.w);

    if sets.w == 0
        sets.comparisonMethod = 'mass_pcc';
    else
        sets.comparisonMethod = 'mpnan';
    end

    % small test. Todo: do not save the whole rezMaxMP since it gives
    % memory issues
%     ticBytes(gcp);
% 
%     [rezMaxMP] = hca_compare_distance(barGen(passingThreshBars),theoryStruct(1:100), sets );
% tocBytes(gcp);

    [rezMaxMP] = hca_compare_distance(barGen(passingThreshBars), theoryStruct, sets );
    


%     cc(idy,idx) = rezMaxMP{1}(idx,idy).maxcoef(1);

%     import Core.export_coefs;
%     export_coefs(theoryStruct,rezMaxMP,bestBarStretchMP,barGen(passingThreshBars),[sets.dirName, '_MP_w=',num2str(sets.w),'_']);
%     save([sets.dirName, num2str(sets.w),'_', num2str(min(sF)),'sf_rez.mat'],'rezMaxMP','passingThreshBars','sets','-v7.3');

    % Get info for all the coefficients. Save this as a separate matrix
    allCoefs = cellfun(@(x) x{1}, rezMaxMP,'un',false);
    matAllCoefs =  cat(3, allCoefs{:});
    save([savedir,'_', num2str(sets.w),'_', num2str(min(sF)),'sf_allcoefs',dateVec,'.mat'],'matAllCoefs','passingThreshBars','sets','-v7.3');
    outputLocs{wIdx}{1} = [savedir, '_',num2str(sets.w),'_', num2str(min(sF)),'sf_allcoefs',dateVec,'.mat'];
    tic
    maxCoef = cell(1,size(matAllCoefs,1));
    maxOr = cell(1,size(matAllCoefs,1));
    maxPos = cell(1,size(matAllCoefs,1));
    maxSecondPos = cell(1,size(matAllCoefs,1));
    maxlen = cell(1,size(matAllCoefs,1));
    bestSF = cell(1,size(matAllCoefs,1));

%% Now find the best coefficient from matAllCoefs (using a cascading or whatever scheme)
    for barid =1:size(matAllCoefs,1)
        [singleCoef , singlePos ] =  max(matAllCoefs(barid,:,:),[],2);
        pos  = squeeze(singlePos)';
        maxCoef{barid} =  squeeze(singleCoef);
        maxOr{barid} = arrayfun(@(x,y) rezMaxMP{x}{3}(barid,y), 1:length(rezMaxMP),pos);
        maxPos{barid} = arrayfun(@(x,y) rezMaxMP{x}{2}(barid,y), 1:length(rezMaxMP),pos);
        maxSecondPos{barid} = arrayfun(@(x,y) rezMaxMP{x}{4}(barid,y), 1:length(rezMaxMP),pos);
        maxlen{barid} = arrayfun(@(x,y) rezMaxMP{x}{5}(barid,y), 1:length(rezMaxMP),pos);
        bestSF{barid} =  sets.theory.stretchFactors(pos);
    end
    toc
  
    matFilepath = export_coefs_local(thryNames,maxCoef,maxOr,maxPos,maxlen, bestSF, barcodeNames(passingThreshBars),[savedir, '_MP_w=',num2str(sets.w),'_'],dateVec);
%     export_coefs(theoryStruct(1:100),rezMax,bestBarStretchMP,barGen(passingThreshBars),[sets.dirName, '_MP_w=',num2str(sets.w),'_']);
    outputLocs{wIdx}{2} = matFilepath;
    outputLocs{wIdx}{3} = sets.thryFile;




end

clear theoryStruct;
% delete(gcp('nocreate'));
    %


t = toc(t0);
% quick_visual_plot(16,9242,barGen,rezMax,bestBarStretch,theoryStruct)

%  super_quick_plot(16,barGen,comparisonStruct,theoryStruct)
% sigmatches = find(allNums ==1)
% for i=1:length(sigmatches)
%     quick_visual_plot(sigmatches(i),9242,barGen,rezMax,bestBarStretch,theoryStruct)
% end

% cell2mat(refNums(sigmatches))
% refNums(signMatch)
% theoryStruct([cell2mat(refNums(signMatch))]).name;

% local_alignment_gui(sets)
end

