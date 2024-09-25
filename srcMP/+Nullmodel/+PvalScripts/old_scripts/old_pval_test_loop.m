% test: generate p-values for bargroups
if ispc
    SCAMP_LINE = '.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else
    sets.SCAMP_LINE = '~/SCAMP/';
    SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end

% delete(gcp('nocreate'))
numWorkers = 30;

% parpool('local',numWorkers)
% This will depend on these parameters: overlap length, number of attempts, re-scaling factors, psf.
% possibly re-scaling factors and psf have minor relevance, then it would
% be mainly the other two parameters that we would need to tune.
%%
% overlap length in general will be fixed for the analysis.
MIN_OVERLAP_PIXELS = 300;


% for the analysis, the total number of fitting attempts will be important,
% not the lengths of individual barcodes. For the analysis we take common
% length, and short barcode
import Nullmodel.gen_random;

NUM_RAND_FRAGMENTS = 100;
PSF_WIDTH_PIXELS = 300/110;
RAND_LENGTH_MIN = 1200;
RAND_LENGTH_2 = 800;

% sF = 1;%

sF = 0.95:0.01:1.05;


% barSynth = cellfun(
% barSynth = cell(1,length(barcodes));
% for i=1:length(barcodes)
%     barSynth{i}.rawBarcode = barcodes{i};
% end

%
%  calculates all MP and MPI
NUM_RAND_FRAGMENTS = 200;
NN = 100;
out ='output';
tic
import Nullmodel.gen_scores;
[maxPCC,totLen] = gen_scores(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,RAND_LENGTH_2,MIN_OVERLAP_PIXELS,NN,sF,out,SCAMP_LINE);
toc

% 
% %% now fit p-val dist. on these. just some tests
% import Zeromodel.beta_ev_fit;
% tic
% [a_fit, n_fit] = beta_ev_fit(maxPCC{1}, [4 1], [inf inf], [4 1], false(1, 2));
% toc
% k=10;
% 
% totLenRescaledCur = totLen(k);
% import Zeromodel.beta_ev_params;
% 
% N_FIXED = totLenRescaledCur*(RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS+1);
% [parameters] = beta_ev_params(maxPCC{k}, 200, N_FIXED);
%     
% [parameters] = beta_ev_params(maxPCC{k}, 200);
% parameters(1)
% % plot fit vs beta ev pdf
% k=5;
% xx=0:0.001:1;
% import Zeromodel.beta_ev_pdf;
% [p] = beta_ev_pdf(xx,  0.1572*MIN_OVERLAP_PIXELS, 1, 2*(500+10*(k-1)+500-2*MIN_OVERLAP_PIXELS));
% figure,plot(xx,p)
% hold on
% histogram(maxPCC{k},'Normalization','pdf')

% now calculate params actually
%      [parameters] = beta_ev_params(maxPCC{k}, MIN_OVERLAP_PIXELS/2);
import Zeromodel.beta_ev_params;

casetest=1;a_fit2=[];n_fit2=[];
for k=1:NN
%     k
%     tic
%     import Zeromodel.beta_ev_params;
    N_FIXED = totLen(k)*(RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS+1);

    if casetest==1
        [parameters] = beta_ev_params(maxPCC{k}, MIN_OVERLAP_PIXELS/2);
    else
        [parameters] = beta_ev_params(maxPCC{k}, 200, N_FIXED);
    end
    %     toc
    a_fit2(k) = parameters(1);
    n_fit2(k) = parameters(2);

%     [a_fit(k), n_fit(k)] = beta_ev_fit(maxPCC{k}, [2 1], [inf inf], [4 (RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS)], [true false]);

end
% a parameters seems to be constant, so we take the mean
aPar = mean(a_fit2);

aPar/MIN_OVERLAP_PIXELS

c1 = polyfit(1:length(a_fit2),a_fit2,1); % fit a line.


% coef=aPar/MIN_OVERLAP_PIXELS

figure,plot(n_fit2)
hold on
kk=1:NN;
plot(2*(RAND_LENGTH_MIN+10*(kk-1)-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_2-MIN_OVERLAP_PIXELS)/300)

% plot(2*(RAND_LENGTH_MIN+10*(kk-1)+RAND_LENGTH_2-2*MIN_OVERLAP_PIXELS))

par2 = (2*(RAND_LENGTH_MIN+10*(kk-1)-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_2-MIN_OVERLAP_PIXELS))./n_fit2
mean(par2)

%% parameter 2.. recalculated based on fixed first parameter
for k=1:NN
    k
%     tic
%     import Zeromodel.beta_ev_params;
%     N_FIXED = totLen(k)*(RAND_LENGTH_MIN+(k-1)*10-MIN_OVERLAP_PIXELS+1);

    m = length(maxPCC{k});
    cc = maxPCC{k};
    x2 = aPar;
    
    denom = (-1 / m * sum(log(1 + betainc(cc.^2, 1/2, x2 / 2 - 1))) + log(2));
    n2(k) = max(0, 1 ./ denom);
end
%%

% now, x1 length is always 500, N2 length os 500+10*(k-1);
% need to remove overlap?
x = 2*(RAND_LENGTH_MIN+10*(kk-1)-MIN_OVERLAP_PIXELS)*(RAND_LENGTH_2-MIN_OVERLAP_PIXELS)/par2;
figure,plot(n_fit2);hold on;plot(x)

% 
% n2test=n2;
% x = 3*(RAND_LENGTH_MIN+([1:100]-1)*10-MIN_OVERLAP_PIXELS+1)+RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS;

% x = 17*(RAND_LENGTH_MIN+([1:100]-1)*10-MIN_OVERLAP_PIXELS+1)+RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS;
% 
% c = polyfit(x,n_fit2,1)
% y_est = polyval(c,x);
% figure,plot(x,n_fit2)
% hold on
% plot(x,y_est,'r--','LineWidth',2)
%
k=1;
xx=0:0.001:1;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx, aPar, 1, x);%y_est(k)
f=figure,plot(xx,p)
hold on
histogram(maxPCC{k},100,'Normalization','pdf')
set(gca,'YScale','log')

intrestingPCC = maxPCC{k}(maxPCC{k}>0.7);
%%
%     N_FIXED = (RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS+1)*(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS+1)*11;
% 
% tic
%  [parameters] = beta_ev_params(maxPCC{1}, 100, N_FIXED);
% parameters(1)
% toc
% % parameters(2)
% tic
%     [a1, n1] = beta_ev_fit(maxPCC{1}, [2 1], [inf inf], [40 N_FIXED], [false true])
% toc
%     [a1, n1] = beta_ev_fit(maxPCC{1}, [2 1], [inf inf], [40 N_FIXED], [false false])

    
% figure,plot([RAND_LENGTH_MIN:10:RAND_LENGTH_MIN+(NN-1)*10]-MIN_OVERLAP_PIXELS+(RAND_LENGTH_MIN-MIN_OVERLAP_PIXELS)*length(sF),n_fit2)
% 
% cellfun(@(x) nanmean(x),maxPCC)
% 
% cell2mat(maxPCC)
% 
% 
% pccMat = cell2mat(maxPCC');


% import TempStuff.gen_random_fragments
% [barcodes, refBarcode] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD);

%% Simpler:
NUM_RAND_FRAGMENTS = 5000;
import Nullmodel.gen_random;
[barcodes] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,0);
barcodeSynth = cell(1,length(barcodes));
for j=1:length(barcodes)
    barcodeSynth{j}.rawBarcode = barcodes{j};
    barcodeSynth{j}.rawBitmask = ones(1,length(barcodeSynth{j}.rawBarcode ));
end
import Nullmodel.gen_random;
[barcodes2] = gen_random(1,PSF_WIDTH_PIXELS,RAND_LENGTH_2,0);


[a,msg] = mkdir(strcat('output',timestamp));
i=1;theoryStruct=[];
theoryStructSynth{i}.filename = fullfile(strcat('output',timestamp),'theory_barcode.txt');
fileID = fopen(theoryStructSynth{i}.filename,'w');
fprintf(fileID,strcat([' %2.' num2str(14) 'f ']), barcodes2{1});
fclose(fileID);
% theoryStruct{i}.meanBpExt_nm = 0.3;
% theoryStruct{i}.psfSigmaWidth_nm = 300;
% theoryStruct{i}.length = length(refBarcode);
% theoryStruct{i}.isLinearTF = 0;
% theoryStruct{i}.name = 'Synthetic theory';

numWorkers = 30;
MIN_OVERLAP_PIXELS = 300;
% minLen = 500; %150kb? or less? depends on application // if GUI, user selects thi
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);
import Core.compare_mp_all
ix = 1;
[~,~,~,~,compStr] = ...
    compare_mp_all(theoryStructSynth,barcodeSynth,minLen,ix, timestamp,sF,MIN_OVERLAP_PIXELS,numWorkers);

%%
pccs= cellfun(@(x) x.maxcoef,compStr);
k=1;
xx=0:0.001:1;
import Zeromodel.beta_ev_pdf;
[p] = beta_ev_pdf(xx, 40, 1, 10000);%y_est(k)
f=figure,%plot(xx,p)
hold on
histogram(pccs,50,'Normalization','pdf')
set(gca,'YScale','log')

%% TEST: pval for real data. Here we take barcodes from unrelated bacteria, and calculate the same.