% 10/01/22 Synthetic
% version of chrom_assembly_script_vs_ref_synthetic
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% STEP1: gen synthetic data

%% gen synthetic
import Nullmodel.gen_random_fragments;

% parameters
TOTAL_RAND_LENGTH = 10000; % total length
PSF_WIDTH_PIXELS = 4.7; % psf
NUM_RAND_FRAGMENTS = 100; % number of random fragments
MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
FRAGMENT_LENGTH_STD = 100; %std of length of fragment
ADDED_NOISE_MEAN = 0.4; % additive noise mean
ADDED_NOISE_STD = 2.0; % additive noise std
FRAGMENT_STRETCH_STD = 0.02; % fragment length-rescale factor std 0.02
IS_CIRCULAR = 1; % whether circular barcode
    
% generate barcodes /also add pcc?
[barcodeGen, refBarcode, origPos, origFlip,origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD,IS_CIRCULAR);


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


% synthStr2 =[];
% for i=1:length(barcodeGen)
%     for j=1:length(barcodeGen)
%         if i~=j
%             pA = synthStr{i}.pos;
%             pB = synthStr{j}.pos;
%             lpA =  synthStr{i}.lengthMatch;
%             lpB =  synthStr{j}.lengthMatch;
%             st = max(pA,pB); % which barcode is more to the left
%             stop = min(pA+lpA-1,pB+lpB-1);
%             synthStr2(i,j).posB = pB-st+1;
%             synthStr2(i,j).posA = pA-st+1; % inaccurate for scaled barcodes..
%             a = [ barcodeGen{i}.rawBarcode  barcodeGen{i}.rawBarcode];
%             b =[ barcodeGen{j}.rawBarcode barcodeGen{j}.rawBarcode];
%             if synthStr{i}.or==1
%                 a = fliplr(a);
%             end
%             if synthStr{j}.or==1
%                 b = fliplr(b);
%             end
%             aFul = a(pA-st+1:pA-st+1+stop-st );
%             bFul =b(-(pB-st+1):-(pB-st+1)+stop-st);
% %             figure,plot(zscore(aFul))
% %             hold on
% %             plot(zscore(bFul))
%       
%             synthStr2(i,j).fulloverlapPosA = pB-st+1:pB+stop-1;
%             synthStr2(i,j).fulloverlapPosARescaled = pA-st+1:pA+stop-1;
%             synthStr2(i,j).fullscore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
% 
%             synthStr2(i,j).overlaplen = length(aFul);
% 
%             synthStr2(i,j).fullscore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
% %             overlapStruct(k,iy).overlaplen = length(aFul);
% %             overlapPos = 
%         else
%              synthStr2(i,j).fullscore =nan;
%         end   
%     end
% end

[a,msg] = mkdir(strcat('output',timestamp));
i=1;theoryStruct=[];
theoryStruct{i}.filename = fullfile(strcat('output',timestamp),'theory_barcode.txt');
fileID = fopen(theoryStruct{i}.filename,'w');
fprintf(fileID,strcat([' %2.' num2str(14) 'f ']), refBarcode);
fclose(fileID);
theoryStruct{i}.meanBpExt_nm = 0.3;
theoryStruct{i}.psfSigmaWidth_nm = 600;
theoryStruct{i}.length = length(refBarcode);
theoryStruct{i}.isLinearTF = 0;
theoryStruct{i}.name = 'Synthetic theory';

sets.comparisonMethod = 'mass_pcc';
sets.filter = 0;
sets.filterSettings.filter = 0;
% nmbp = barcodeGen{1}.nmbp;
nmpx = 110; % 208?
psf = 300; %300 nm


bpPx = 500;
import CBT.Hca.UI.Helper.plot_best_pos;
fig1 = plot_best_pos([], synthStr, [], [], [],theoryStruct{1}.length,1);%1/bpPx*10^6
grid on
grid minor

%

% %%
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% 
% %  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
% %      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)
% 
% import Thry.gen_theoretical;
% [theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);

%% calc PCC main calc // alternative mp.
sF = 0.95:0.01:1.05;
minOverlap = 150;
tic
import Core.calc_overlap_pcc_sort_m;
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap);
toc

% alternative: MP
tic
% [overlapStructMP] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap);
[overlapStructMP] = calc_overlap_mp(barcodeGen',sF, minOverlap,timestamp)
toc

%% now p-values
import Zeromodel.pval_for_struct;
pscores = pval_for_struct(overlapStruct,0.6,0.05);

%% Similarity test: calculate std for N best match. This will be used to filter out
% bars locally matching to several.
import Nullmodel.local_sim_test;
[curStd, pairs] = local_sim_test(pscores, barcodeGen', overlapStruct,minOverlap,150,3000,timestamp,overlapStructMP);
%%%%

%% create graph
 psc = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1

    filtM = triu(psc);
    filtM(filtM==0) = nan;


import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(pscores,sum(~isnan(filtM(:))),barcodeGen,curStd,100,overlapStruct);

%% coverage map
import Core.create_graph_coverage_map;
[data,prev] = create_graph_coverage_map(finalgraph,overlapStruct,3)

%% visual
PCC_OVERLAP = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1

filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');
idxPair = 6;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,xId, yId,barStruct);
%%
import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,1,2,barStruct);


save('synthdatatest.mat','curStd','data','prev','finalgraph','Ggraphs','pscores','overlapStruct','barcodeGen','sF','minOverlap','refBarcode','synthStr', '-v7');