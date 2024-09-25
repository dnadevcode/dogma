% TEST to compare

rootdir =  '/media/albyback/4EF8DB03F8DAE86B/DATA/POSTDOC/sub/';

import Core.create_barcodegen;
[barcodeGen, kymoStructs, sets] = create_barcodegen(rootdir,1,10,1,timestamp,1);


minLen = 350; %150kb? or less? depends on application // if GUI, user selects this
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(barLens>minLen);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false);...
    cellfun(@(x) x.nmbp,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask','nmbp'},2);

nmbp = 0.192;
% nmbp = 0.245;
nmpx = 110;
% mean(unique(cellfun(@(x) x.nmbp, barcodeGen)))
% nmbp = 0.22;

%%
fastaFile =  '/home/albyback/git/bargroupingprototype/test/DA32087.fasta';
%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)

import Thry.gen_theoretical;
[theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);


sF = 0.8:0.01:1.2;
numWorkers = 2;

% parpool('local',numWorkers)
% addpath(genpath('C:\Users\Lenovo\git\bargroupingprototype'))
% cd C:\Users\Lenovo\git\bargroupingprototype

sets.SCAMP_LINE = '/home/albyback/postdocData/test_transloc/SCAMP/';

if ispc % for PC get from the initialization
    SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else %
    SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end

MIN_OVERLAP_PIXELS = 250;

foldSynth = 'barcoli';

[namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);

% names2{1} = ;
names2 = arrayfun(@(x) theoryStruct{1}.filename, 1:length(namesBar),'un',false);
tic
import Core.calc_overlap;
[mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,names2,namesBar);
toc

%% Check the re-scaling factors along the barcode. 
strFac = cell(1,length(mp1));
pos = cell(1,length(mp1));
bestSTF = zeros(1,length(mp1));

for i=1:length(mp1)
    [a,b] = sort(mp1{i},'Desc');
    strFac{i} = stridx{i}(b(1:(sum(stridx{i}==1)-MIN_OVERLAP_PIXELS )));
    pos{i} = mpI1{i}(b(1:(sum(stridx{i}==1)-MIN_OVERLAP_PIXELS )));
    bestSTF(i)= sF(abs(strFac{i}(1)));
end
    
% check if all position along the genome have the same re-scaling factor. 
% MP can be extended to full overlap, and we find best re-scaling factor
% for all those positions
% The same procedure can be done for the all to all comparison.    


%% for experiment
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
%%
% we compare xId with yId. yId is un-rescaled and xId is re-scaled
% xId = 195
% yId = 163
% we find the positions on mp1 and mpI1
positionsBar = find(baridx2{xId}==yId);

mpBar = mp1{xId}(positionsBar);
mpIBar = mpI1{xId}(positionsBar);
rescaleBar = mpI1{xId}(positionsBar);

mpIBar = mpIBar(mpIBar>-1);


% figure2:
strFactors = stridx{xId}(mpIBar+1);


% figure 1: mpI positions
% figure,plot(mpIBar)
% 
% import Core.get_two_bars_alignment_params;
% [bar1Name, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{xId}, mpI1{xId}, ps, baridx2{xId},sF);

% should only keep the re-scaling factors for the pixels where the barcodes
% are overlapping
 sortedidx2 = find(pksUnique1{xId}==yId);
%     sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);

pos = LOCS1{xId}(pksUniquePos1{xId}(sortedidx2)); % position on MP.
import Core.get_two_bars_alignment_params;

[~, pA, pB,rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{xId},mpI1{xId},pos, baridx2{xId},sF); 
% pos1 = -pA+1:-pA+length(bRescaled);
% pos2 = -pB+1:-pB+length(b2);


%

vec = [yId xId];
% [bar1idx,bar2idx] = [29],[49];
bar1idx=vec(2);
bar2idx=vec(1);
sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);
import Core.plot_best_match;
[f,pos1, bar1, pos2, bar2,pcc1,pcc2] = plot_best_match(sortedidx2,stridx,mpI1{bar1idx},LOCS1{bar1idx},pksUniquePos1{bar1idx},bar1idx,baridx2{bar1idx},sF,barStruct,MIN_OVERLAP_PIXELS);
pcc1
pcc2

[C,IA,IB] = intersect(pos1,pos2);
xlim([C(1) C(end)])

nexttile([1 2]),plot(-pB:-pB-1+length(strFactors),[strFactors])
xlim([C(1) C(end)])
nexttile([1 2]),plot(-pB:-pB-1+length(strFactors),[mpBar(~isnan(mpBar))])
xlim([C(1) C(end)])

% xlim([-pA+1 -pA+length(bar1)])

% todo: estimate PSF/NMBP/NMPX based on data


%% re-scaling compare test newer structure