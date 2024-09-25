% using MP, we locally compare bars


fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta','/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta','/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta'};

%  [t,matFilepathShort,theoryStruct, sets,theoryGen] = ...
%      HCA_theory_parallel(theory_names,meanBpExt_nm,setsFile,theoryfile,theoryfold)

psffac = 1;
nmbp = 0.225;
nmpx = 110*psffac;
import Thry.gen_theoretical;
[theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmpx);


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

for idxTheory = 1:2;
k=1;
numWorkers = 30;
out = timestamp;
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
figure,plot(MPthr{1})
% figure,plot(MPthr{2})

[a,b] = max(MPthr{1});
loc = MPIthr1{1}(b)+1;

bar1 = thr{1}(b:b+MIN_OVERLAP_PIXELS-1);
bar2 = thr{1}(loc:loc+MIN_OVERLAP_PIXELS-1);

figure,plot(bar1)
hold on
plot(bar2)
zscore(bar1,1)*zscore(bar2,1)'/length(bar1)

b1 = zscore(bar1,1);
b2 =zscore(bar2,1);
% 1-sum((b1-b2).^2)/(2*length(bar1))
% 1-sum((b1-b2).^2)/(2*length(bar1))

%% both

theoryStruct2{1}.filename = 'flip.txt';
writematrix(fliplr(importdata( theoryStruct{1}.filename )),theoryStruct2{1}.filename );

% idxTheory = 2;
k=1;
MIN_OVERLAP_PIXELS = 150;
numWorkers = 30;
out = timestamp;
if ispc
    command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=.\'...
       theoryStruct2{1}.filename  ' --input_a_file_name=.\'...
       theoryStruct{2}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
       ' --output_a_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_mp' ])) ...
       ' --output_a_index_file_name=.\' fullfile(out,strcat([num2str(k) '_' num2str(k) '_index' ])) ]);
else
     
       command = strcat([SCAMP_LINE ' --window=' num2str(MIN_OVERLAP_PIXELS) ' --input_a_file_name=\'...
        theoryStruct2{1}.filename ' --input_a_file_name=\'...
        theoryStruct{2}.filename ' --num_cpu_workers=' num2str(numWorkers) ' --no_gpu --output_pearson --print_debug_info'...
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

figure,plot(mp1{1})
