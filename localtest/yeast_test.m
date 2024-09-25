load('/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/barcodeGen993.mat')
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out)

barcodeGen1  = barcodeGen;

minLen = 200;
barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen1);
barcodeGen1 = barcodeGen1(barLens>=minLen);
lengths1 = cellfun(@(x) sum(x.rawBitmask),barcodeGen1);


[barcodeGenNonCirc,nonper] = filt_periodic(barcodeGen1,timestamp,1,3);%0.8:0.01:1.2
% barcodeGen = barcodeGenNonCirc;

minOverlap = 200;
sF = 1;
tic
import Core.calc_overlap_pcc_sort;
[overlapStruct] = calc_overlap_pcc_sort([barcodeGenNonCirc], sF,minOverlap);

%%
%% check dist for non matching real
% overlapStruct1 = overlapStruct(end/2+1:end,1:end/2);
overlapStruct1 = overlapStruct;


% rel = triu(overlapStruct1);

allpairs = overlapStruct1(:);

rgg = 200:10:450;
ppvec =[];
for numPts = rgg;
    sc = [];
    rPos = randi(numel(overlapStruct1),1,10000);

    for i=1:length(rPos)
        if length(allpairs(rPos(i)).xcorrs)~=1
            newpts =allpairs(rPos(i)).xcorrs(allpairs(rPos(i)).numElts{1}==numPts);
            if length(newpts(:))==4
                sc =[ sc; newpts(1) ];
            end
        end
    end
    import Zeromodel.beta_ev_nu;
    [pars] = beta_ev_nu(sc, numPts/3);
    
    ppvec = [ppvec pars/numPts];
    pars/numPts
end

% figure,plot(rgg,ppvec.*(300:10:800))
% figure,plot(rgg,ppvec.*rgg./(rgg+200))
figure,plot(rgg,ppvec)
ylabel('\nu/overlap length')
xlabel('Overlap pixels')
cnu = polyfit(rgg,ppvec,1);

predPSF =1/(mean(ppvec)*sqrt(2*pi))

%%

PCC_OVERLAP = reshape([overlapStruct1.sc], size(overlapStruct1,1),size(overlapStruct1,2));
overlaplen = reshape([overlapStruct1.overlaplen], size(overlapStruct1,1),size(overlapStruct1,2));
bestBarStretch = reshape([overlapStruct1.bestBarStretch], size(overlapStruct1,1),size(overlapStruct1,2));

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGenNonCirc,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGenNonCirc,'un',false)]',{'rawBarcode','rawBitmask'},2);

idxPair = 4;
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
% normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(filtM(:),'desc','MissingPlacement','last');

[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));

% length(sortedVals(sortedVals>sortedVals2(2)))
% number of false positives based on two dataset, can then compare with
% pvals


import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct1,xId, yId,barStruct);
%

% sort out bars with repetitive pattern
%%




fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/sequence1.fasta',...
    '/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/sequence2.fasta',...
    '/proj/snic2022-5-384/users/x_albdv/data/Yeast/new/sequence3.fasta'};

% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};

% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)

barcodeGen{1}.nmbp = 0.25;
psffac = 1;
nmbp = barcodeGen{1}.nmbp*psffac;
nmpx = 254*psffac;
import Thry.gen_theoretical;
[theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);

sets.comparisonMethod = 'mass_pcc';

sF = 0.9:0.01:1.1;
[comparisonStruct] = compare_to_t(barcodeGen,theoryStruct,sF,sets)


bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], comparisonStruct, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

% figure,plot(compStr{1}.allposA-compStr{1}.allposB)

theoryStruct2{1} = theoryStruct{1};

theoryStruct2{1}.filename = 'combined_barcode.txt';
theoryStruct2{1}.name = 'combined_barcode';
longBar = [];
for i=1:length(theoryStruct)
    barcode = importdata(theoryStruct{i}.filename );
    longBar = [longBar barcode nan];
end

fd = fopen(theoryStruct2{1}.filename,'w' );
fprintf(fd,'%5.3f ',longBar);
fclose(fd);


[compStr] = compare_to_t_mp(barcodeGen,theoryStruct2,1,150)

% idx=5;
N = 50;
pos = [];
stdPos =[];
for idx=1:length(compStr)
    pos{idx} = (compStr{idx}.allposA-compStr{idx}.allposB);
    stdPos(idx) = std(pos{idx}(1:min(end,N)));
end
bG = barcodeGen(stdPos<5);

% figure,plot(pos)

mean(cellfun(@(x) x.maxcoef(1),comparisonStruct))

%%

% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta','/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)
% 
% psffac = 1;
% nmbp = 0.225;
% nmpx = 110*psffac;
% import Thry.gen_theoretical;
% [theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);


%%

if ispc % for PC get from the initialization
    sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP' ;
    SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else %
%     sets.SCAMP_LINE = '/home/albyback/postdocData/test_transloc/SCAMP/';
    sets.SCAMP_LINE = '~/SCAMP/';

    SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end
%%
MIN_OVERLAP_PIXELS = 150;
k=1;
numWorkers = 30;
out = timestamp;
mkdir(out)
for idxTheory = 1:3;

%
if ispc
    command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
       theoryStruct{idxTheory}.filename  ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
else
     
       command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
        theoryStruct{idxTheory}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
end
[~,message] = system(command); %--output_pearson

   MPIthr1{idxTheory} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));
    
    %also get PCC scores.
%     tic
    % a bit slower import to make sure imports correctly on win
fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
raw2 = textscan(fid, '%s ');
fclose(fid);
nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
MPthr{idxTheory} = nan(length(nonanValues),1);
MPthr{idxTheory}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
thr{idxTheory} = importdata(   theoryStruct{idxTheory}.filename);

end
%     toc
% figure,plot(mp1{1})
idx=3;

figure,
tiledlayout(2,1)
nexttile
plot(MPthr{idx})
nexttile
plot(MPIthr1{idx})

% figure,plot(MPthr{2})

% idx=2;
[a,b] = max(MPthr{idx});
loc = MPIthr1{idx}(b)+1;

bar1 = thr{idx}(b:b+MIN_OVERLAP_PIXELS-1);
bar2 = thr{idx}(loc:loc+MIN_OVERLAP_PIXELS-1);

figure,plot(bar1)
hold on
plot(bar2)
zscore(bar1,1)*zscore(bar2,1)'/length(bar1)

b1 = zscore(bar1,1);
b2 =zscore(bar2,1);
% 1-sum((b1-b2).^2)/(2*length(bar1))
% 1-sum((b1-b2).^2)/(2*length(bar1))

%% both
ix1=1;
ix2=2;

% idxTheory = 2;
k=1;
MIN_OVERLAP_PIXELS = 150;
numWorkers = 30;
out = timestamp;
if ispc
    command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
       theoryStruct{ix1}.filename  ' --input_a_file_name=.\'...
       theoryStruct{ix2}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
else
     
       command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
        theoryStruct{ix1}.filename ' --input_a_file_name=\'...
        theoryStruct{ix2}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
end
[~,message] = system(command);

   mpI1{k} = importdata(fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])));
    
    %also get PCC scores.
%     tic
    % a bit slower import to make sure imports correctly on win
fid = fopen(fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])));
raw2 = textscan(fid, '%s ');
fclose(fid);
nonanValues = cellfun(@(x) x(1)~='-',raw2{1});
mp1{k} = nan(length(nonanValues),1);
mp1{k}(nonanValues) = sscanf(sprintf(' %s',raw2{1}{nonanValues}),'%f');
%    

idx = 1;
figure,
tiledlayout(2,1)
nexttile
plot(mp1{idx})
nexttile
plot(mpI1{idx})