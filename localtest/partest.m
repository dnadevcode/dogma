c = parcluster('tetralith');
c.AdditionalProperties.AccountName = 'snic2022-5-384';
c.AdditionalProperties.WallTime = '02:00:00';
c.AdditionalProperties.Reservation = '';
c.AdditionalProperties.NumNodes = 1;
c.AdditionalProperties.ProcsPerNode = 1; 

c.saveProfile

% https://www.hpc.iastate.edu/guides/using-matlab-parallel-server
tic
myjob = batch(c,'mywave','AutoAddClientPath',false)
wait(myjob)
% diary(myjob)
toc

c.AdditionalProperties.ProcsPerNode = 9; 
tic
myjob = batch(c,'parallel_mywave','pool', 8, 'AutoAddClientPath',false);
wait(myjob)
% diary(myjob)
toc


c.AdditionalProperties.NumNodes = 1;
c.AdditionalProperties.ProcsPerNode = 16; 
c.AdditionalProperties.Reservation = '';
tic
myjob = batch(c,'parallel_eigen','pool', c.AdditionalProperties.NumNodes*(c.AdditionalProperties.ProcsPerNode-1), 'AutoAddClientPath',false);
wait(myjob)
% 
diary(myjob)
toc

%%

c.AdditionalProperties.NumNodes = 1;
c.AdditionalProperties.ProcsPerNode = 30; 
c.AdditionalProperties.Reservation = '';


testSet ={ '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/',...
    '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_Old E.coli_ EF365_after shrink finder/'};
numF = 10;
testSet2 = {'/proj/snic2022-5-384/users/x_albdv/data/Yeast/','/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/'};

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

minLen = 250;
synth = 1;

import Zeromodel.prep_data;
% [barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

% tic
myjob = batch(c,'Zeromodel.prep_data',4,{testSet2,minLen,minLen,synth},'pool', c.AdditionalProperties.NumNodes*(c.AdditionalProperties.ProcsPerNode-1), 'AutoAddClientPath',false,'AttachedFiles',{'hca_parallel_settings.txt','default_edge_detection.txt'});
wait(myjob)
[barcodeGen1,barcodeGen2,lengths1,lengths2] = myjob.fetchOutputs{:};


minOverlap = 200;
sF = 0.8:0.01:1.2;
c.AdditionalProperties.NumNodes = 2;
c.AdditionalProperties.ProcsPerNode = 32; 

myjob3 = batch(c,'Core.calc_overlap_pcc',1,{[barcodeGen1 barcodeGen2], sF,minOverlap},'pool', c.AdditionalProperties.NumNodes*(c.AdditionalProperties.ProcsPerNode-1), 'AutoAddClientPath',false);
wait(myjob3)

[overlapStruct] = myjob3.fetchOutputs{:};
% tic
% import Core.calc_overlap_pcc;
% [overlapStruct] = calc_overlap_pcc([barcodeGen1 barcodeGen2], sF,minOverlap);
% toc

load(myjob,'barcodeGen1','barcodeGen2','len1','len2');

%% simple overlap
minLen2 = 700;
minLen1 = 10000;

import Zeromodel.prep_data;
% [barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

% tic
myjob = batch(c,'Zeromodel.prep_data',4,{testSet2,minLen1,minLen2,synth},'pool', c.AdditionalProperties.NumNodes*(c.AdditionalProperties.ProcsPerNode-1), 'AutoAddClientPath',false,'AttachedFiles',{'hca_parallel_settings.txt','default_edge_detection.txt'});
wait(myjob)
[barcodeGen1,barcodeGen2,lengths1,lengths2] = myjob.fetchOutputs{:};


minOverlap = 300;
sF = 0.8:0.01:1.2;
c.AdditionalProperties.NumNodes = 2;
c.AdditionalProperties.ProcsPerNode = 32; 

myjob3 = batch(c,'Core.calc_overlap_pcc',1,{[barcodeGen1 barcodeGen2], sF,minOverlap},'pool', c.AdditionalProperties.NumNodes*(c.AdditionalProperties.ProcsPerNode-1), 'AutoAddClientPath',false);
wait(myjob3)
[overlapStruct] = myjob3.fetchOutputs{:};

%%

synth=1;
c.AdditionalProperties.NumNodes = 1;
c.AdditionalProperties.ProcsPerNode = 30; 
minLen2 = 700;
minLen1 = 700;

import Zeromodel.prep_data;
% [barcodeGen1,barcodeGen2,lengths1,lengths2] = prep_data(testSet,numF,minLen,synth)

% tic
myjob = batch(c,'Zeromodel.prep_data',4,{testSet2,minLen1,minLen2,synth},'pool', c.AdditionalProperties.NumNodes*(c.AdditionalProperties.ProcsPerNode-1), 'AutoAddClientPath',false,'AttachedFiles',{'hca_parallel_settings.txt','default_edge_detection.txt'});
wait(myjob)
[barcodeGen1,barcodeGen2,lengths1,lengths2] = myjob.fetchOutputs{:};

c.AdditionalProperties.NumNodes = 4;
c.AdditionalProperties.ProcsPerNode = 32; 

minOverlap = 300;
sF = 1;

myjob5 = batch(c,'Core.calc_overlap_pcc_null',1,{[barcodeGen1 barcodeGen2], sF,minOverlap},'pool', c.AdditionalProperties.NumNodes*(c.AdditionalProperties.ProcsPerNode-1), 'AutoAddClientPath',false);

[overlapStruct] = myjob5.fetchOutputs{:};

% tic
% import Core.calc_overlap_pcc_null;
% [overlapStruct] = calc_overlap_pcc_null([barcodeGen1 barcodeGen2], sF,minOverlap);
% 
