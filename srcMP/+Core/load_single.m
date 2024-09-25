function [barcodeGen,passingKymos,kymoStructs,kS,filelist] = load_single(foldL)

rootdir = foldL;
filelist = dir(fullfile(rootdir, '**/*.tif'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

% arrayfun(@(x) ~isempty(strfind('Lambda',filelist(x).folder)),1:length(filelist),'un',true)

sets.kymosets.kymofilefold = arrayfun(@(x) filelist(x).folder,1:length(filelist),'un',false);
sets.kymosets.filenames = arrayfun(@(x) filelist(x).name,1:length(filelist),'un',false);
sets.output.matDirpath = 'output';
% single fold
% tiffs = dir(fullfile(foldL,'*.tif'));

sets.edgeDetectionSettings.method = 'Zscore';
sets.alignMethod = 0;
sets.filterSettings.filter = 0;
sets.skipEdgeDetection = 0;
sets.bitmasking.untrustedPx = 5;
sets.genConsensus = 0;

import Core.load_barcodes_avg;
[barcodeGen,passingKymos,kymoStructs,kS] = load_barcodes_avg(sets,0,1);

end

