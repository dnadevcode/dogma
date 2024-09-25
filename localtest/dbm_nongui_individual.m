
useGUI = 0;

import OldDBM.General.SettingsWrapper;
defaultSettingsFilepath = SettingsWrapper.get_default_newDBM_ini_filepath();
if not(exist(defaultSettingsFilepath, 'file'))
    defaultSettingsFilepath = '';
end
dbmOSW = SettingsWrapper.import_dbm_settings_from_ini(defaultSettingsFilepath);

dbmOSW.DBMSettingsstruct.dbmtool = 'hpfl-odm'; 
dbmOSW.DBMSettingsstruct.askForDBMtoolSettings = 0;

dbmOSW.DBMSettingsstruct.movies.askForMovies = 0;

dbmOSW.DBMSettingsstruct.detectlambdas = 0;
dbmOSW.DBMSettingsstruct.initialAngle = 0; % used if moleculeAngleValidation = 0
dbmOSW.DBMSettingsstruct.maxLambdaLen = inf;
dbmOSW.DBMSettingsstruct.moleculeAngleValidation = 1;

% files = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220318\tif\*.tif');
% files = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220112\tif\*.tif');

Folder   = 'C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\';
FileList = dir(fullfile(Folder, '**', '*.czi'))

% tiffname = 'Experiment-5868';
tiffname = 'Experiment-5692';

idx = find(arrayfun(@(x) ~isempty(strfind(FileList(x).name,tiffname)),1:length(FileList)));

% files = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220201\tif\*.tif');
dbmOSW.DBMSettingsstruct.outputDirpath = 'C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\test\';

% files = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\ECOLIMOV\*.tif');
% dbmOSW.DBMSettingsstruct.outputDirpath = 'C:\Users\Lenovo\postdoc\DATA\Chromosome\ECOLIMOV\kymos\';


filesC = arrayfun(@(x) fullfile(FileList(x).folder,FileList(x).name),idx,'un',false);
filesC{1} = strrep(filesC{1},'\czi\','\tif\');
filesC{1} = strrep(filesC{1},'.czi','.tif');

% filesC = arrayfun(@(x) fullfile(files(x).folder,files(x).name),15:20,'un',false);

dbmOSW.DBMSettingsstruct.movies.movieNames = filesC;

dbmOSW.DBMSettingsstruct.max_number_of_frames = 10;
dbmOSW.DBMSettingsstruct.max_f = 10;
dbmOSW.DBMSettingsstruct.minLen = 200;

% dbmOSW.DBMSettingsstruct.movies.movieNames = {filesC{1:3}};

% fd =fopen(dbmOSW.DBMSettingsstruct.movies.movieFile);
% filePh = fopen(dbmOSW.DBMSettingsstruct.movies.movieFile,'w');
% fprintf(filePh,'%s\n',filesC{:});
% fclose(filePh);
 
dbmOSW.DBMSettingsstruct.genome_assembly_pipeline = 1;
dbmOSW.DBMSettingsstruct.auto_run = 0;

useGUI = 0;

dbmOSW.DBMSettingsstruct.choose_output_folder = 0;
% DBM_Gui(useGUI,dbmOSW)
dna_barcode_matchmaker(0,dbmOSW);
 
%% czi to tif
% data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220318\czi\*5835*.czi');
% % % % 
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220111\czi\*.czi');
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220112\czi\*.czi');
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220127\czi\*.czi');
DBM4.convert_czi_to_tif(data,0);

data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220318\czi\*.czi');
% 
test.bargrouping_convert_czi_to_tif(data,0);
%
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220221\czi\*.czi');
% 
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220127\czi\*.czi');
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220127\czi\*.czi');
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220208\czi\*.czi');
test.bargrouping_convert_czi_to_tif(data,0);
%  
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220203\czi\*.czi');
test.bargrouping_convert_czi_to_tif(data,0);
%  
data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\Yeast_Luis_processed\220201\czi\*.czi');
test.bargrouping_convert_czi_to_tif(data,0);

% data = dir("C:\Users\Lenovo\postdoc\DATA\Calibration\fluorsegmen_project\2019-12-13 experiments\2019-12-13 lungcancercells\DAPI\FOV2_DAPI_gainVariation\nd2\*.nd2");
% % % 
% DBM4.convert_czi_to_tif(data,0);

data = dir('C:\Users\Lenovo\postdoc\DATA\Chromosome\czi files\czi files\*.czi');
% test.bargrouping_convert_czi_to_tif(data,0);
DBM4.convert_czi_to_tif_direct(data,0);
