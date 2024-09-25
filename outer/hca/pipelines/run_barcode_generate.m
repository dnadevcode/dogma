% load default settings
function barcodeGen=run_barcode_generate(tiffname,minLen)

if nargin < 1
    tiffname = 'C:\Users\Lenovo\git\lldev\OutputFiles\RawKymos\yeast_Data\*.tif';
    sets.minLen = 50; % miniimum length
    sets.maxLen = 200;
end


% load default settings
import CBT.Hca.Import.import_hca_settings;
[sets] = import_hca_settings('hca_settings.txt');
sets.kymosets.askforkymos = 0;
sets.kymosets.askforsets = 0;


% load tiff names // could be directly laoded to 
tiffs = dir(tiffname);

sets.kymosets.filenames = arrayfun(@(x) x.name, tiffs, 'UniformOutput', false);
sets.kymosets.kymofilefold = arrayfun(@(x) x.folder, tiffs, 'UniformOutput', false);

% load kymo names into settings
% import CBT.Hca.Settings.get_user_settings;
% sets = get_user_settings(sets);

% timestamp for the results
sets.timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
sets.alignMethod = 3;


% add kymographs
import CBT.Hca.Import.add_kymographs_fun;
[kymoStructs] = add_kymographs_fun(sets);

%  put the kymographs into the structure
import CBT.Hca.Core.edit_kymographs_fun;
kymoStructs = edit_kymographs_fun(kymoStructs,sets.timeFramesNr);

% align kymos
import CBT.Hca.Core.align_kymos;
[kymoStructs] = align_kymos(sets,kymoStructs);

% generate barcodes
import CBT.Hca.Core.gen_barcodes;
barcodeGen =  CBT.Hca.Core.gen_barcodes(kymoStructs, sets);

