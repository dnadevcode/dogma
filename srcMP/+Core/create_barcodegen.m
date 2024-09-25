function [barcodeGen,kymoStructs, sets] = create_barcodegen(fold,type,timeFramesNr,alignMethod,timestamp,timeframest)

%   Args: 
%       fold - folder with the data
%       type - whether lambda or kymograph    
%       timeFramesNr - number of timeframes 
%       alignMethod - alingment method

%       timeframest - starting timeframe

% fold = '/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/all';

% timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

if nargin < 5
    timeframest = 1;
end
%% test 1) compare individual kymographs against theory.
import CBT.Hca.Import.import_hca_settings;
[sets] = import_hca_settings('hca_parallel_settings.txt');

sets.kymosets.kymoFile = 'kymos.txt';

% check if there are .mat files with barcode
filelist = dir(fullfile(fold, '**/*.mat'));  %get list of files and folders in any subfolder
if ~isempty(filelist)
    struc = load(fullfile(filelist(1).folder,filelist(1).name));
    try % tru to see if exists such a field
        barcodeGen = struc.barcodeGen;
        kymoStructs = [];
        sets = [];
        sets.timestamp = timestamp;
        mkdir(timestamp);
        
        return;
    catch
    end
end

filelist = dir(fullfile(fold, '**/*kymograph.tif'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list
% 


name = 'lambda';
% keep specific kind of kymographs
lambdaF = arrayfun(@(x) isempty(strfind(lower(filelist(x).folder),name)),1:length(filelist));
notLambdaFolds = find(lambdaF);
lambdaFolds =  find(~lambdaF);
% [strfind('kymographs',filelist.folder)]

try
    nmbpidx = arrayfun(@(x) strfind(lower(filelist(x).folder),'nmperbp'),1:length(filelist));
    nmbp = arrayfun(@(x) str2num(filelist(x).folder(nmbpidx(x)-5:nmbpidx(x)-1)) ,1:length(filelist));
    nmbpnotLambdaFolds = nmbp(find(lambdaF));
catch
    nmbpnotLambdaFolds = 0.2*ones(1,length(notLambdaFolds));
end
    

fd = fopen(sets.kymosets.kymoFile,'w');
if type==1
    arrayfun(@(x) fprintf(fd,'%s\n',fullfile(filelist(x).folder,filelist(x).name)),notLambdaFolds);
else
    arrayfun(@(x) fprintf(fd,'%s\n',fullfile(filelist(x).folder,filelist(x).name)),lambdaFolds);
end
fclose(fd);

sets.timestamp = timestamp;

import CBT.Hca.Settings.get_user_settings;
sets = get_user_settings(sets);

%%
sets.timeFramesNr = timeFramesNr;
sets.alignMethod = alignMethod;
%mkdir(sets.timestamp);
        

%% add kymos
matKymopathShort = fullfile(sets.output.matDirpath, strcat(['kymos_' sprintf('%s_%s', timestamp) '.txt']));
fd = fopen(matKymopathShort,'w');    fclose(fd);
                  
% predefine structure
kymoStructs = cell(1,length(sets.kymosets.filenames));
    
for i=1:length(kymoStructs)
    kymoStructs{i}.name = sets.kymosets.filenames{i};
    kymoStructs{i}.unalignedKymo = flipud(imread(fullfile(sets.kymosets.kymofilefold{i},kymoStructs{i}.name)));  
    try
        kymoStructs{i}.unalignedBitmask = flipud(imread(fullfile(sets.kymosets.kymofilefold{i},strrep(kymoStructs{i}.name,'kymograph','bitmask')))); 
    catch
    end
    fd = fopen(matKymopathShort,'a'); fprintf(fd, '%s \n',fullfile(sets.kymosets.kymofilefold{i},kymoStructs{i}.name)); fclose(fd);
end
     
% import CBT.Hca.Import.add_kymographs_fun;
% [kymoStructs] = add_kymographs_fun(sets);

for i=1:length(kymoStructs)
    kymoStructs{i}.unalignedKymo = kymoStructs{i}.unalignedKymo(timeframest:sets.timeFramesNr+timeframest-1,:);
    try
        kymoStructs{i}.unalignedBitmask = kymoStructs{i}.unalignedBitmask(timeframest:sets.timeFramesNr+timeframest-1,:);
    catch
    end
end 
        
%         
% %  put the kymographs into the structure
% import CBT.Hca.Core.edit_kymographs_fun;
% kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr,timeframest);

% align kymos
import CBT.Hca.Core.align_kymos;
[kymoStructs] = align_kymos(sets,kymoStructs);

% generate barcodes
import CBT.Hca.Core.gen_barcodes;
barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);

for i=1:length(barcodeGen)
    barcodeGen{i}.nmbp = nmbpnotLambdaFolds(i);
end

