% 10/03/22 real (yeast)
% version of chrom_assembly_script_vs_ref_synthetic
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% STEP1: gen exp data
caseM = 1;

% yeast
if caseM == 1
    load('/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/barcodeGen993.mat');
    kymonames = cellfun(@(x) x.name,barcodeGen,'un',false)';
%     find(cellfun(@(x) ~isempty(strfind(x,'5589')),kymonames))
    load('/proj/snic2022-5-384/users/x_albdv/data/Yeast/test/barcodeGenTest.mat')

else
    testSet = {'/proj/snic2022-5-384/users/x_albdv/data/CHR/Mapping_Old E.coli_ EF365_after shrink finder/'};
testSet = {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/Mapping_New E.coli/testA/'};

testSet = {'C:\Users\Lenovo\postdoc\DATA\Chromosome\ECOLIMOV\kymos\'};

    import Zeromodel.prep_data;
    [barcodeGen,barcodeGen2,lengths1,lengths2] = prep_data(testSet,10,500,0)
end


out=strcat('output',timestamp);
mkdir(out)
% barcodeGen  = barcodeGen;


% only keep barcodes longer than minimum length.
minLen = 200; %150kb? or less? depends on application
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
barcodeGen = barcodeGen(barLens>minLen);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);



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


%% calc PCC main calc // alternative mp.
sF = 0.85:.01:1.15;
% sF = 1;
minOverlap = 150;
tic
import Core.calc_overlap_pcc_sort_m;
[overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap);

% % alternative: MP
tic
[overlapStructMP] = calc_overlap_mp(barcodeGen,sF, minOverlap,timestamp)
toc
%% now p-values/ should convert to Stouffer Z-score
import Zeromodel.pval_for_struct;
pscores = pval_for_struct(overlapStruct,0.6,0.03);
NN = sum(~isnan(pscores));

%% Similarity test: calculate std for N best match. This will be used to filter out
% bars locally matching to several.
import Nullmodel.local_sim_test;
[curStd,pairs,scrs] = local_sim_test(pscores, barcodeGen, overlapStruct,minOverlap,150,3000,timestamp,overlapStructMP);
%%%%
% MIN_OVERLAP_PIXELS,N,NN,timestamp,overlapStructMP)
%% create graph
import Core.create_overlap_graph
[finalgraph,Ggraphs] = create_overlap_graph(pscores,351,barcodeGen,curStd,100,overlapStruct);

%% coverage map
import Core.create_graph_coverage_map;
[data,prev] = create_graph_coverage_map(finalgraph,overlapStruct,1)

%% visual
PCC_OVERLAP = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1

PCC_OVERLAP = reshape([overlapStruct.sc], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1
% PCC_OVERLAP = reshape([overlapStructMP.score], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1

% pcc2=PCC_OVERLAP(2:end,1:end-1);

% for i=1:size(PCC_OVERLAP,2)-1
%     PCC_OVERLAP(i,i+1) = nan;
% %     PCC_OVERLAP(i-1,i) = nan;
% end

filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');
%%
idxPair = 25;
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,xId, yId,barStruct);

% import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, overlapStructMP,xId, yId);
%%
save('realdatatest1.mat','curStd','data','prev','finalgraph','Ggraphs','pscores','overlapStruct','barcodeGen','sF','minOverlap', '-v7');
% save('realdatatest2.mat','curStd','data','prev','finalgraph','Ggraphs','pscores','overlapStruct','barcodeGen','sF','minOverlap','refBarcode','synthStr', '-v7');