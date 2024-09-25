% 09/06/22 True positive test

% version of chrom_assembly_script_vs_ref_synthetic
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% STEP1: gen synthetic data

%% gen synthetic
import Nullmodel.gen_random_fragments;

% parameters
TOTAL_RAND_LENGTH = 10000; % total length
PSF_WIDTH_PIXELS = 2.7; % psf
NUM_RAND_FRAGMENTS = 300; % number of random fragments
MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
FRAGMENT_LENGTH_STD = 100; %std of length of fragment
ADDED_NOISE_MEAN = 0.4; % additive noise mean
ADDED_NOISE_STD = 1.0; % additive noise std
FRAGMENT_STRETCH_STD = 0.02; % fragment length-rescale factor std 0.02
IS_CIRCULAR = 1; % whether circular barcode
    
% generate barcodes /also add pcc?
[barcodeGen,refBarcode,origPos, origFlip,origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD,IS_CIRCULAR);


% only keep barcodes longer than minimum length.
minLen = 300; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(barLens>minLen);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);


% plot the fragment positions. Requires hca function
synthStr = cell(1,length(barcodeGen));
for i=1:length(barcodeGen)
    synthStr{i}.idx = 1;
    synthStr{i}.pos = origPos(i);
    synthStr{i}.lengthMatch = round(length(barcodeGen{i}.rawBarcode)*origStr(i));
    synthStr{i}.or = origFlip(i);
    synthStr{i}.rf = origStr(i);
end

[a,msg] = mkdir(strcat('output',timestamp));
i=1;theoryStruct=[];
theoryStruct{i}.filename = fullfile(strcat('output',timestamp),'theory_barcode.txt');
fileID = fopen(theoryStruct{i}.filename,'w');
fprintf(fileID,strcat([' %2.' num2str(14) 'f ']), refBarcode);
fclose(fileID);
theoryStruct{i}.meanBpExt_nm = 0.3;
theoryStruct{i}.psfSigmaWidth_nm = 300;
theoryStruct{i}.length = length(refBarcode);
theoryStruct{i}.isLinearTF = 0;
theoryStruct{i}.name = 'Synthetic theory';

sets.comparisonMethod = 'mass_pcc';
sets.filter = 0;
sets.filterSettings.filter = 0;
% nmbp = barcodeGen{1}.nmbp;
nmpx = 110; % 208?
psf = 300; %300 nm


% %%
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% 
% %  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
% %      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)
% 
% import Thry.gen_theoretical;
% [theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);

%% calc PCC main calc
sF = 0.8:0.01:1.2;
minOverlap = 150;
tic
import Core.calc_overlap_pcc_sort_m;
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap);

%% PCC
sF = 0.95:0.01:1.05;

rezMax=[];bestBarStretch=[];bestLength=[];
for i=1:length(theoryStruct)
    tic
    import CBT.Hca.Core.Comparison.on_compare;
    [rezMax{i},bestBarStretch{i},bestLength{i}] = on_compare(barcodeGen,...
        theoryStruct{i},sets.comparisonMethod, sF,[],50,[],sets.filterSettings);
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
% import CBT.Hca.UI.get_display_results;
% [res] = get_display_results(barcodeGen,[], comparisonStruct, theoryStruct, sets);
%% MP
numWorkers = 30;
MIN_OVERLAP_PIXELS_BL = 300;
minLen_BL = 300;
% minLen = 500; %150kb? or less? depends on application // if GUI, user selects thi
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen');
[a,b] = find(barLens>minLen_BL);
% find(b==44)
import Core.compare_mp_all
ix = 1;
[mpI1individual,mp1individual,maxMP,stridxindividual,compStr] = ...
    compare_mp_all(theoryStruct,barcodeGen',minLen_BL,ix, timestamp,sF,MIN_OVERLAP_PIXELS_BL,numWorkers);

lengthBorders = cumsum(cellfun(@(x) x.length,theoryStruct));
% lengthBorders = lengthBorders(ix);
% % fig1 = figure;
% bpPx = 500;
% import CBT.Hca.UI.Helper.plot_best_pos;
% fig1 = plot_best_pos([], compStr, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
% grid on
% grid minor
%     saveas(fig1,'figs/fig4.png')
% 
% % 
% plot_best_pos([], comparisonStruct, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

%% Overlaps
MIN_OVERLAP_PIXELS=300;
barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);
% have to save separately..
foldSynth= 'barcoli';
% have to save separately..
[namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);
% skip one barcode
[names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);



if ispc % for PC get from the initialization
    sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP' ;
    SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else %
%     sets.SCAMP_LINE = '/home/albyback/postdocData/test_transloc/SCAMP/';
    sets.SCAMP_LINE = '~/SCAMP/';

    SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end


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
% import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, overlapStruct,11,6);

PCC_OVERLAP = reshape([overlapStruct.score], size(overlapStruct,1),size(overlapStruct,2));
overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));
% lenOverlap = 

%%
%% filter scores based on best overlap
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');

% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 1;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));


% PCC_OVERLAP(xId,yId)

import Core.plot_match_simple;
[f] = plot_match_simple(barStruct, overlapStruct,xId,yId);

overlapStruct(xId,yId).pB-overlapStruct(xId,yId).pA+1

% [res] = get_display_results(barcodeGen([xId yId]),[], comparisonStruct([xId yId]), theoryStruct, sets);
%   import CBT.Hca.UI.Helper.plot_best_pos;
  bpPx= 500;
   f= plot_best_pos([], compStr([find(b==xId) find(b==yId)]), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%        saveas(f,'figs/fig5.png')
% 
   plot_best_pos([], comparisonStruct([xId yId]), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
% %     plot_best_pos([], comparisonStruct, [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
% 
% 
% pthresh = 0.00001;
% import Core.simple_bargroup;
% [orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,yId,pthresh,1);
% % [orBars,refac,pA,pB,bars] =simple_bargroup(barStruct, overlapStruct,7,pthresh,1);
% 
%  f=   plot_best_pos([], comparisonStruct(bars), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);
%         saveas(f,'figs/fig7.png')
%   
%  plot_best_pos([], compStr(bars), [], [], [],theoryStruct{1}.length,1/bpPx*10^6);

%%
nT = 1000;%length(sortedIds);
tp = zeros(1,nT);
fp = zeros(1,nT);
tp2 = zeros(1,nT);
fp2 = zeros(1,nT);
gt =  zeros(1,nT);
for ii=1:nT
    idxPair = ii;
    [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    px1 = comparisonStruct{xId}.pos(1):comparisonStruct{xId}.pos(1)+comparisonStruct{xId}.lengthMatch-1;
    px2 = comparisonStruct{yId}.pos(1):comparisonStruct{yId}.pos(1)+comparisonStruct{yId}.lengthMatch-1;
    overLap = intersect(mod(px1,theoryStruct{1}.length),mod(px2,theoryStruct{1}.length));
    if length(overLap)>=100
        tp(ii)=1;
    else
        fp(ii)=1;
    end

        px1 = compStr{xId}.pB(1)-compStr{xId}.pA(1)+1:compStr{xId}.pB(1)-compStr{xId}.pA(1)+1+compStr{xId}.lengthMatch-1;
        px2 = compStr{yId}.pB(1)-compStr{yId}.pA(1)+1:compStr{yId}.pB(1)-compStr{yId}.pA(1)+1+compStr{yId}.lengthMatch-1;
    overLap = intersect(mod(px1,theoryStruct{1}.length),mod(px2,theoryStruct{1}.length));
    if length(overLap)>=100
        tp2(ii)=1;
    else
        fp2(ii)=1;
    end

        px1 = synthStr{xId}.pos(1):synthStr{xId}.pos(1)+synthStr{xId}.lengthMatch-1;
    px2 = synthStr{yId}.pos(1):synthStr{yId}.pos(1)+synthStr{yId}.lengthMatch-1;
    overLap = intersect(mod(px1,theoryStruct{1}.length),mod(px2,theoryStruct{1}.length));
    if length(overLap)>=100
%         tp(ii)=1;
    else
        gt(ii)=1;
    end

end

f=figure,plot(cumsum(fp2))
hold on
plot(cumsum(fp))
plot(cumsum(gt))

% Unrecognized function or variable 'tp1'.

legend({'MP','PCC','GT'})
saveas(f,'figs/fig8.png')
