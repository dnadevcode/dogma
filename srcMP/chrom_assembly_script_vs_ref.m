%% whether to ask user for inputs. 1- yes, 0 no 
GUI = 1;
dataFold = 'C:\Users\Lenovo\postdoc\DATA\'; % kymo folder
fastaFold = pwd; % fasta folder

%%

% add hca to path for theory calculation and pcc comparison
if ispc
    addpath(genpath('C:\Users\Lenovo\git\hca'));
end
% addpath(genpath('/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/hca'))
addpath(genpath(pwd));

% create a time-stamp for collecting the results
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

% command-line root-dir for tifs
if ispc
    rootdir = 'C:\Users\Lenovo\postdoc\DATA\Frag\all2';
    rootdir = 'C:\Users\Lenovo\postdoc\DATA\Alma-20220405T082558Z-001\Alma\3 BACs AFF1 gene';
    rootdir = 'C:\Users\Lenovo\postdoc\DATA\Alma-20220405T082558Z-001\Alma\3 BACs KMT2A gene';
    rootdir='C:\Users\Lenovo\postdoc\DATA\Alma-20220405T082558Z-001\Alma\3 BACs MLLT3 gene';
    
%     rootdir = 'C:\Users\Lenovo\postdoc\DATA\Alma-20220405T082558Z-001\Alma\3 BACs KMT2A gene\Raw kymographs 1 BAC';
%     rootdir = 'C:\Users\Lenovo\postdoc\DATA\Alma-20220405T082558Z-001\Alma\3 BACs KMT2A gene\Raw kymographs 2 BACs';

else
    rootdir = '/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/all';
end

% possible root dirs for different examples
if GUI
    rootdir = uigetdir(dataFold);
end


%% test 1) compare individual kymographs against theory.

import Core.create_barcodegen;
[barcodeGen, kymoStructs, sets] = create_barcodegen(rootdir,1,10,1,timestamp);
disp(strcat(['Detected ' num2str(length(barcodeGen)) ' barcodes']));

% only keep barcodes longer than minimum length.
minLen = 200; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(barLens>minLen);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);

% create theory..
if GUI
    [fastaFile,pathname] = uigetfile({'*.fasta','*.*'},'MultiSelect','on',dataFold); limits = [];
    fastaFile = fullfile(pathname,fastaFile);
else
    % fastaFile =  'DA32087.fasta';
    fastaFile = '018_final_polish.fasta'; limits = [];
    fastaFile = 'C:\Users\Lenovo\postdoc\DATA\seq\chr11.fna'; gap = 100000; limits = [87856154-gap 88062206+gap];  
    fastaFile = 'C:\Users\Lenovo\postdoc\DATA\seq\chr10.fna'; gap = 100000; limits = [21524616-gap 21743630+gap];
    %KMT2A
    fastaFile = 'C:\Users\Lenovo\postdoc\DATA\seq\chr11.fna'; gap = 500000; limits = [118436492-gap 118526832+gap];%KMT2A
    %MLTT3
    fastaFile = 'C:\Users\Lenovo\postdoc\DATA\seq\chr9.fna'; gap = 500000; limits = [20341669-gap 20622499+gap];%KMT2A
end

% pixelWidth_nm,psfSigmaWidth_nm,isLinearTF,resultsDir
nmbp = 0.205; 
nmbp = 0.193;
nmbp = 0.22;
nmpx = 130;
nmpsf = 300;
% n
    sets.comparisonMethod = 'mass_pcc';

    sets.filter = 0;
if GUI
    import Gui.get_cam_pars;
    [nmbp, nmpx, nmpsf] = get_cam_pars();
    sets.filter = 1;
    sets.filterMethod = 1;
    sets.filterSize = 50;
end

%%
% nmbp = 0.26;
% nmbp,pixelWidth_nm,psfSigmaWidth_nm,isLinearTF,resultsDir
[theoryStruct] = gen_theory_for_chrom(nmbp,nmpx,nmpsf,1, sets.timestamp,limits, fastaFile);
% [theoryStruct] = gen_theory_for_chrom(nmbp,'018_final_polish.fasta');

% figure,plot(importdata(theoryStruct{1}.filename))

% RUN comparison
% todo: compare in the same way as MP. 
% can create one long file from theoryStruct.
rezMax=[];bestBarStretch=[];bestLength=[];
for i=1:length(theoryStruct)
 tic
 import CBT.Hca.Core.Comparison.on_compare;
[rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGen,...
    theoryStruct{i},sets.comparisonMethod,0.9:0.01:1.1,[],50,[],sets);
toc
end

comparisonStructAll = rezMax;
for i=1:length(comparisonStructAll)
    for j=1:length(bestBarStretch{i})
        comparisonStructAll{i}{j}.bestBarStretch = bestBarStretch{i}(j);
        comparisonStructAll{i}{j}.length = bestLength{i}(j);
    end
end
import CBT.Hca.Core.Comparison.combine_theory_results;
[comparisonStruct] = combine_theory_results(theoryStruct, rezMax,bestBarStretch,bestLength);

sets.timeFramesNr = nan;
sets.displayResults=1
sets.userDefinedSeqCushion = 0;
sets.genConsensus = 0;
import CBT.Hca.UI.get_display_results;
[res,fig1] = get_display_results(barcodeGen,[], comparisonStruct, theoryStruct, sets);
% saveas(fig1,'transloc_fig1.eps','epsc')
% create consensus
[outConsensus,coverage, pval,consensus,f] = gen_reference_based_assembly(barcodeGen,comparisonStruct,theoryStruct,timestamp);
% saveas(f,'transloc_fig2.eps','epsc')

% 
% save('MLLT3_consensus.mat','consensus','limits','gap','barcodeGen','comparisonStruct','theoryStruct','timestamp')

% save('KMT2A_consensus.mat','consensus','limits','gap','barcodeGen','comparisonStruct','theoryStruct','timestamp')
% figure,plot(consensus)
% hold on
% plot(outConsensus(1,:))
% posGene = [gap/(nmpx/nmbp) gap/(nmpx/nmbp)+(diff(limits)-2*gap)/(nmpx/nmbp)];
% plot(posGene,ones(1,length(posGene)),'blackx-')
%% consensus specific maps
%
barIdx = [8 30];
[outConsensus] = gen_reference_based_assembly(barcodeGen(barIdx),comparisonStruct(barIdx),theoryStruct);

outConsensus2=circshift(outConsensus,[0,10]);
f=figure,plot(find(~isnan(nanmean(outConsensus2(2:end,:)))),outConsensus2(2:end,~isnan(nanmean(outConsensus2(2:end,:))))'); 
