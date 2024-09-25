addpath(genpath('/home/albyback/git/bargroupingprototype'))


GUI = 0; % add GUI elements
tFrames = 1;
% foldTest = 'C:\Users\Lenovo\postdoc\DATA\Mapping_New E.coli\Mapping_New E.coli\final data\20220610_Sample DA32087-4-st1_110nmPERpx_0.192nmPERbp\kymos';
% foldTest = 'C:\Users\Lenovo\postdoc\DATA\Mapping_New E.coli\Mapping_New E.coli\final data\20220613_Sample DA32087-4-st2_110nmPERpx_0.193nmPERbp\kymos';
% foldTest = '/media/albyback/4EF8DB03F8DAE86B/DATA/POSTDOC/E. coli Assembly - data for Lund/Mapping_New E.coli/final data/20220613_Sample DA32087-4-st2_110nmPERpx_0.193nmPERbp/kymos/';
foldTest = '/media/albyback/4EF8DB03F8DAE86B/DATA/POSTDOC/E. coli Assembly - data for Lund/Mapping_New E.coli/final data/20220610_Sample DA32087-4-st1_110nmPERpx_0.192nmPERbp/kymos/';

foldTest = '/media/albyback/4EF8DB03F8DAE86B/DATA/POSTDOC/Final_Data_Shrink_Sorted/';
foldTest = '/media/albyback/4EF8DB03F8DAE86B/DATA/POSTDOC/sub/';

foldTest = '/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/';
if ispc % add hca to path for theory calculation and pcc comparison
    foldTest = 'C:\Users\Lenovo\postdoc\DATA\sub\Mapping_New E.coli\Final_Data_Shrink_Sorted/';
end
    % SCAMP_LINE = 
%%
if ispc % add hca to path for theory calculation and pcc comparison
    addpath(genpath('C:\Users\Lenovo\git\hca'))
end
% addpath(genpath('/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/hca'))

% addpath(genpath(pwd));
% dataFold = pwd;
% create a time-stamp for collecting the results

% if ~exist('bargroupingprototype','dir')
addpath(genpath(pwd));
dataFold = pwd;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
mkdir(strcat('output',timestamp));
save(fullfile('pipelines','currentTimestampbg.mat'),'timestamp');
% else
%     load('currentTimestampbg.mat')
% end
       
% timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
% mkdir(strcat('output',timestamp));

%% test 1) compare individual kymographs against theory.
import CBT.Hca.Import.import_hca_settings;
[sets] = import_hca_settings('hca_parallel_settings.txt');

sets.kymosets.kymoFile = 'kymos.txt';

% first load barcodes
% fold = 'C:\Users\Lenovo\postdoc\DATA\Frag\Fragment för matchning - Data till Erik\Mapping_New E.coli\final data';
if ispc
    fold='C:\Users\Lenovo\postdoc\DATA\Frag\all';
else
    fold = '/pfs/proj/nobackup/fs/projnb10/snic2021-5-132/all';
end

if GUI
    fold = uigetdir(dataFold);
end

try
    fold = foldTest;
catch 
end
import Core.create_barcodegen;
[barcodeGen, kymoStructs, sets] = create_barcodegen(fold,1,tFrames,1,timestamp,1);
disp(strcat(['Detected ' num2str(length(barcodeGen)) ' barcodes']));

% choose barcode gen good

minLen = 350; %150kb? or less? depends on application // if GUI, user selects this
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(barLens>minLen);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);

%% Comparison against theory, beginning of the script is identical, but then follows _vs_ref
% 
%% ASSEMBLY MP
% delete(gcp('nocreate'))
% numWorkers = 4;
% 
% parpool('local',numWorkers)
% addpath(genpath('C:\Users\Lenovo\git\bargroupingprototype'))
% cd C:\Users\Lenovo\git\bargroupingprototype

% sets.SCAMP_LINE = '/home/albyback/postdocData/test_transloc/SCAMP/';

if ispc % for PC get from the initialization
    sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP\';
    SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else %
%     SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
    sets.SCAMP_LINE = '~/SCAMP/' ;
        SCAMP_LINE = strcat([sets.SCAMP_LINE 'build/SCAMP']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows


end

MIN_OVERLAP_PIXELS = 250;
% NMPX = 110;
% NMPSF = 200;
% NMBP = 0.22;
sF = 0.95:0.01:1.05; %0.85:0.01:1.15;
% sF = 0.85:0.01:1.15;

%% save barcodex into txts

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
% have to save separately..
foldSynth= 'barcoli';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);


% skip one barcode
[names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);


%%

%%  calculates all MP and MPI
numWorkers = 30;
tic
import Core.calc_overlap;
[mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,names2,namesBar);
toc


for ii=1:length(baridx2)
    baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
end
% we want to create a nicer structure for all-to-all comparison and
% contains easily accesible data.
tic
import Core.mp_res_to_struct;
[overlapStruct] = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
toc
import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, overlapStruct,11,6);
%%

NN=length(mpI1);
PKS1=cell(1,NN);
LOCS1=cell(1,NN);
pksUnique1 = cell(1,NN);
pksUniquePos1=cell(1,NN);
barcodePair1=cell(1,NN);
rescaleFactorPair1=cell(1,NN);
orPair1=cell(1,NN);
pvals =zeros(1,NN);
for k=1:NN
%     k
    [PKS1{k},LOCS1{k},pksUnique1{k},pksUniquePos1{k},barcodePair1{k},rescaleFactorPair1{k},orPair1{k}] = ...
    Core.unique_peak_barcodes(mp1{k},mpI1{k},baridx2,stridx,k);
end
%%
% calculate full overlap PCC scores
import Core.mp_to_full_overlap;
[PCC_OVERLAP, len1, len2,lenOverlap,PCC_MP] = mp_to_full_overlap(stridx,pksUnique1,pksUniquePos1,LOCS1,mpI1,baridx2, MIN_OVERLAP_PIXELS, barStruct,sF);

%%
import Core.overlap_graph
% G = overlap_graph
[pval,G] = overlap_graph(MIN_OVERLAP_PIXELS,PCC_OVERLAP,len1,len2)


import Core.layout_consensus
layout_consensus

%% single bargroup.
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:).*sqrt(lenOverlap(:));
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');

% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 2;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

PCC_OVERLAP(xId,yId)
bar1idx=xId;bar2idx=yId;
sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);
import Core.plot_best_match;
[f,pos1, bar1, pos2, bar2,pcc1,pcc2] = plot_best_match(sortedidx2,stridx,mpI1{bar1idx},LOCS1{bar1idx},pksUniquePos1{bar1idx},bar1idx,baridx2{bar1idx},sF,barStruct,MIN_OVERLAP_PIXELS);
pcc1
pcc2
%%
% xId = 104
pthresh = 0.00001;
import Core.single_bargroup;
 [barMat,bars,orBars,reFac,pval] = single_bargroup(xId,pthresh,PCC_OVERLAP,PCC_MP,len1,len2,lenOverlap,...
        stridx,mpI1,LOCS1,pksUnique1,pksUniquePos1,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,1)
    
    
% bargroup map
import Core.bargroup_map;
 [barMat,shiftVec,barM,barMzscored] = bargroup_map(bars, reFac, orBars,barcodeGen,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers);
% statistics of the chance of getting random match

% coverage for this bargroup.
figure,plot(sum(~isnan(barM(:,1:end-1))))

    %% overlalp graph bargroup
     pthresh = 0.0001;
    import Core.bargroup_assembly;
     [Ggraphs,numComp,numNodes] = bargroup_assembly(PCC_OVERLAP,PCC_MP, lenOverlap, len1,len2, pksUnique1,sortedidx2,stridx,mpI1,LOCS1,...
    pksUniquePos1,bar1idx,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS,barcodeGen,pthresh);

G = Ggraphs{end};
figure;
plot(G,'Layout','force','ArrowSize',10,'MarkerSize',1)
% 
    import Core.bargroup_consensus;
    
    
    % todo: plot path between two barcodes based on G
    import Core.plot_path

%%
[a,b] = sort(maxMP,'descend');
%%
% import Core.plot_best_match;
ix1 = b(5)
% % ix1= 3
% ix1=195;
ix = 1;
pos = LOCS1{ix1}(pksUniquePos1{ix1}(ix)); % position on MP.

import Core.plot_best_match;
[f,~,~,~,~,~,bar1Name] = plot_best_match(ix,stridx,mpI1{ix1},LOCS1{ix1},pksUniquePos1{ix1},ix1,baridx2{ix1},sF,barStruct,MIN_OVERLAP_PIXELS);

% example_stretch
% example stretch
%%
saveas(f,'fig1.png')
barIdx = [ix1  bar1Name];

barIdx = [1 9];
% if theory generated..
[outConsensus] = gen_reference_based_assembly(barcodeGen(barIdx),comparisonStruct(barIdx),theoryStruct);

outConsensus2=circshift(outConsensus,[0,10]);
f=figure,plot(find(~isnan(nanmean(outConsensus2(2:end,:)))),outConsensus2(2:end,~isnan(nanmean(outConsensus2(2:end,:))))'); 

%% plot specific % 13,59, 29, 49
% vec =[13 59];
vec = [29 49];
% vec = [29 13];
vec = [xId,yId];

vec = [166 225];
% [bar1idx,bar2idx] = [29],[49];
bar1idx=vec(1);
bar2idx=vec(2);
sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);
import Core.plot_best_match;
[f,pos1, bar1, pos2, bar2,pcc1,pcc2] = plot_best_match(sortedidx2,stridx,mpI1{bar1idx},LOCS1{bar1idx},pksUniquePos1{bar1idx},bar1idx,baridx2{bar1idx},sF,barStruct,MIN_OVERLAP_PIXELS);
pcc1
pcc2
saveas(f,'localvsglobal3.fig')
saveas(f,'localvsglobal3.png')

% pcc
pccval = mp1{bar1idx}(LOCS1{bar1idx}(pksUniquePos1{bar1idx}(sortedidx2)))
% pval
% %     % will also compute p-value
a = 0.122*MIN_OVERLAP_PIXELS;
len2 = length(bar1);%lengths(barcodePair1{k}(1));
len1 = length(bar2);
n = 0.004*(len2-MIN_OVERLAP_PIXELS+1)*(len1-MIN_OVERLAP_PIXELS+1);
% 
pval = 1-Zeromodel.beta_ev_cdf(pccval, a, 1, n, 1)



    
    
 %%
allbars = b;
[c,idPos] = sort(b);
allb = ones(1,length(allbars));
%
bvec= [];

graphMatrix = zeros(length(barcodeGen),length(barcodeGen));

for ii=1:10
    if sum(allb) ==0
        break;
    end
    % now loop through first few
    btemp = allbars(find(allb));
    idx = 1;
    bar1 = btemp(idx);

    % take PCC vals for unique matching barcodes
    allPCCvals = PKS1{bar1}(pksUniquePos1{bar1}); %% CHECK/POSSIBLE BUG
    
    % convert to PVAL
    len2 = lengths(pksUnique1{bar1});
    len1 = lengths(bar1)
    pvalPar1 = 0.13*MIN_OVERLAP_PIXELS;
    pvalPar2 = 2*(len1-MIN_OVERLAP_PIXELS)*(len2-MIN_OVERLAP_PIXELS)/100;
    import Zeromodel.beta_ev_pdf;
%     pval =

%     [p] = beta_ev_pdf(allPCCvals', pvalPar1, 1, pvalPar2);%y_est(k)
    for j=1:length(allPCCvals)
        [p(j)] =  1-Zeromodel.beta_ev_cdf(allPCCvals(j), pvalPar1, 1, pvalPar2(j), 1);
    end

    
    % find number above thresh
%     numGoodMatch = allPCCvals > 0.8; % here for PCC
    numGoodMatch = p < 0.01;
    % threshold peaks based on intensity 
    % % now for a single barcode check the cluster barcodes
    % bar1= goodIdx(8);  % 709
    % bar1= 548;  % 709
    peaksToTry = pksUnique1{bar1}(numGoodMatch);
    
    graphMatrix(bar1,peaksToTry) = 1; % just the graph matrix
%     graphMatrix(peaksToTry,bar1) = 1; % just the graph matrix

%     allb(peaksToTry)
    % now create bargroup for bar1.
%     import Core.plot_bargroup;
%     [barMat{ii},bars{ii}] = plot_bargroup(bar1,stridx,mpI1(peaksToTry),LOCS1(peaksToTry),pksUnique1(peaksToTry), pksUniquePos1(peaksToTry),peaksToTry,baridx2,sF,barStruct,MIN_OVERLAP_PIXELS);
     %                                       ix,stridx,mpI                   ,LOCS                   ,pksUnique1                  ,pksUniquePos                   ,k,               baridx, sF, barcodeGenGood, h

    allb(idPos([bar1 peaksToTry])) = 0;
%     bvec{ii} = [bar1 peaksToTry];
end
    figure,imagesc(graphMatrix)
