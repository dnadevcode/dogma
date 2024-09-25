function [kymoStructs,barGen,nmpx,nmbp] = load_kymo_data(dirName,depth,ix,iy,sF,timeFramesNr)

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

% smart extract nmPerbp from folder name
[nmbp, nmpx] = Input.extract_extension_params(sets.dirName);

% spltName = strsplit(sets.dirName ,'_');
% spltName2 = strsplit(spltName{end},'nm');
% spltName3 = strsplit(spltName{end-1},'nm');
% sets.nmbp = str2double(spltName2{1});
% nmpx = str2double(spltName3{1});


% load theory
% thryFileIdx = find(arrayfun(@(x) ~isempty(strfind(thryFiles(x).name,num2str(nmpx))),1:length(thryFiles)));
% sets.thryFile = fullfile(thryFiles(thryFileIdx).folder,thryFiles(thryFileIdx).name);

%%

files = dir(fullfile(sets.dirName,'kymos','*.tif'));

sets.kymosets.filenames = arrayfun(@(x) files(x).name,1:length(files),'un',false);
sets.kymosets.kymofilefold = arrayfun(@(x) files(x).folder,1:length(files),'un',false);

% simdata = 0;
sets.output.matDirpath = 'output';
sets.filterSettings.filter = 0;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 6;
sets.minLen = 150;
sets.genConsensus  = 0;


if nargin < 6
%  following "Strain-level bacterial typing directly from patient
% samples using optical DNA mapping"
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
import Core.load_kymo_data;
[kymoStructs,barGen] = load_kymo_data(sets);

% rescale
import Core.rescale_barcode_data;
[barGen] = rescale_barcode_data(barGen,sets.theory.stretchFactors);


end

