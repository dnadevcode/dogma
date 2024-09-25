


% % test: generate p-values for bargroups
% if ispc
%     SCAMP_LINE = '.\..\SCAMP\build\Release\SCAMP.exe'; % windows
% else
%     sets.SCAMP_LINE = '~/SCAMP/';
%     SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
% end

% delete(gcp('nocreate'))
% numWorkers = 32;
% parpool('local',numWorkers)
%%
NUM_RAND_FRAGMENTS = 500;
PSF_WIDTH_PIXELS = 300/110;
RAND_LENGTH_MIN = 1500;
RAND_LENGTH_2 = 800;

% sF = 1;%

sF = 0.8:0.01:1.2;

import Nullmodel.gen_random;
[~,barcodeGen] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,0);
%

minOverlap = 300;
tic
import Core.calc_overlap_pcc;
[overlapStruct] = calc_overlap_pcc(barcodeGen, sF,minOverlap);
toc

PCC_OVERLAP = reshape([overlapStruct.score], size(overlapStruct,1),size(overlapStruct,2));
overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));
bestBarStretch = reshape([overlapStruct.bestBarStretch], size(overlapStruct,1),size(overlapStruct,2));

% lenOverlap = 

import Zeromodel.beta_ev_params;

% fit the parameters:
maxpcc = PCC_OVERLAP(~isnan(PCC_OVERLAP));
[parameters] = beta_ev_params(maxpcc, minOverlap);

parameters(1)/minOverlap
parameters(2)
% 2*(RAND_LENGTH_MIN+RAND_LENGTH_MIN-2*minOverlap)

%%
overlap = overlaplen(~isnan(PCC_OVERLAP));
f1 = figure,histogram(overlap(:),RAND_LENGTH_MIN-minOverlap+1)
saveas(f1,'fig2.eps','epsc')

% pd = fitdist(r,'exponential');
f=figure;
h = histfit(overlap(:)-minOverlap+1,50,'exponential')
xlabel('Overlap length')
saveas(f,'fig2.eps','epsc')

parameters1 =[];
parameters2 =[];

for j=minOverlap:RAND_LENGTH_MIN
    scores = maxpcc(overlap==j);
    if length(scores) > 100
        import Zeromodel.beta_ev_params;

        [parameters] = beta_ev_params(scores, minOverlap/3);
        parameters1 = [parameters1 parameters(1)];
        parameters2 = [parameters2 parameters(2)];
    end
end
x =minOverlap:minOverlap+length(parameters1)-1;

%     f = figure,plot(x,parameters1)
% y
    c = polyfit(x,parameters1,1)
y_est = polyval(c,x);
f=figure,plot(x,parameters1)
hold on
plot(x,y_est,'r--','LineWidth',2)
saveas(f,'fig4.eps','epsc')

    c = polyfit(x,parameters2,1)
y_est = polyval(c,x);
f=figure,plot(x,parameters2)
hold on
plot(x,y_est,'r--','LineWidth',2)
saveas(f,'fig5.eps','epsc')
    % fit the parameters:
%     maxpcc = PCC_OVERLAP(~isnan(PCC_OVERLAP));


%         h = histfit(scores,30,'gev')

% f=figure;
% h = histfit(scores,50,'gev')
% saveas(f,'fig3.eps','epsc')

% end
% h(1).FaceColor = [.8 .8 1];
% h(2).Color = [.2 .2 .2];
% set(gca,'YScale','log')
% f2 = figure,histogram(PCC_OVERLAP(:))
% % saveas(f2,'h1.png')

%%
 barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
    % have to save separately..
    [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);
    % skip one barcode
    [names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);


    % delete(gcp('nocreate'))
    numWorkers = 32;
    tic
    import Core.calc_overlap;
    [mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,names2,namesBar);
    toc
    
    
    for ii=1:length(baridx2)
    baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
    end
    % we want to create a nicer structure for all-to-all comparison and
    % contains easily accesible data.
%     tic
    import Core.mp_res_to_struct;
    overlapStruct2 = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF,barStruct);
%     toc


%% filter scores based on best overlap
filtM = triu(PCC_OVERLAP);
filtM(filtM==0) = nan;
normalizedScore = filtM(:);%.*sqrt(overlaplen(:)); % score adjusted by sqrt of length

% normalizedScore = filtM(:).*sqrt(overlaplen(:)); % score adjusted by sqrt of length
[sortedVals, sortedIds] = sort(normalizedScore,'desc','MissingPlacement','last');
% 


barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

% [maxScore,maxLoc] = max(normalizedScore);
idxPair = 1
[xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));




import Core.plot_match_pcc;
[f] = plot_match_pcc(barStruct, overlapStruct,xId, yId);
saveas(f,'fig1.eps','epsc')


% % we have
% lenA = RAND_LENGTH_MIN;
% minOverlap = 300;
% 
% % parpool('local',numWorkers)
% % This will depend on these parameters: overlap length, number of attempts, re-scaling factors, psf.
% % possibly re-scaling factors and psf have minor relevance, then it would
% % be mainly the other two parameters that we would need to tune.
% %%
% % overlap length in general will be fixed for the analysis.
% MIN_OVERLAP_PIXELS = 300;
% 
% 
% % for the analysis, the total number of fitting attempts will be important,
% % not the lengths of individual barcodes. For the analysis we take common
% % length, and short barcode
% import Nullmodel.gen_random;
% 
% 
% 
% 
% % barSynth = cellfun(
% % barSynth = cell(1,length(barcodes));
% % for i=1:length(barcodes)
% %     barSynth{i}.rawBarcode = barcodes{i};
% % end
% 
% %
% %  calculates all MP and MPI
% NUM_RAND_FRAGMENTS = 200;
% NN = 100;
% out ='output';
% tic
% import Nullmodel.gen_scores;
% [maxPCC,totLen] = gen_scores(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,RAND_LENGTH_2,MIN_OVERLAP_PIXELS,NN,sF,out,SCAMP_LINE);
% toc
% 

%%
%%
k=1;
par1 = []; % parameter values
par2 = [];

RAND_LENGTH_MIN = 1200;
RAND_LENGTH_2 = 800;

rangePSF = 80:10:170;


% sF = 1;%
% sF = 0.95:0.01:1.05;
rangeSF = [];
for i=0:25
    rangeSF{i+1}=1-i/100:0.01:1+i/100;
end

NUM_RAND_FRAGMENTS = 60;
minOverlap = 300; % always fixed

rangeLen = 350:10:1200;
import Nullmodel.gen_random;
% k=1;
rLen = 1200;
nmpx=110;
sF = {0.8:0.01:1.2};
%% 
% for nmpx=rangePSF % 3 loops each for parameter
% for sF=rangeSF % 3 loops each for parameter
for rLen = rangeLen

    PSF_WIDTH_PIXELS = 300/nmpx;
    RAND_LENGTH_MIN= rLen;
    [~,barcodeGen] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,0);
    %
    
    tic
    import Core.calc_overlap_pcc;
    [overlapStruct] = calc_overlap_pcc(barcodeGen, sF{1},minOverlap);
    toc
    
    PCC_OVERLAP = reshape([overlapStruct.score], size(overlapStruct,1),size(overlapStruct,2));
    overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));
    bestBarStretch = reshape([overlapStruct.bestBarStretch], size(overlapStruct,1),size(overlapStruct,2));
    
    % lenOverlap = 
    
    import Zeromodel.beta_ev_params;
    
    % fit the parameters:
    maxpcc = PCC_OVERLAP(~isnan(PCC_OVERLAP));
    [parameters] = beta_ev_params(maxpcc, minOverlap/3);

    par1(k) = parameters(1)/minOverlap
    par2(k) = parameters(2)
    k=k+1
% 2*(RAND_LENGTH_MIN+RAND_LENGTH_MIN-2*minOverlap)
end
%%
f=figure,
tiledlayout(2,1)
nexttile
plot(rangePSF,par1)
ylabel('\nu')
nexttile
plot(rangePSF,par2)
ylabel('\lambda')
xlabel('nm/px')
saveas(f,'pcc1.png')
save('pcc1.mat','rangePSF','par1','par2')

figure,
tiledlayout(2,1)
nexttile
plot(par1)
ylabel('\nu')
nexttile
plot(par2)
ylabel('\lambda')
xlabel('sF')
saveas(f,'pcc2.png')
save('pcc2.mat','rangeSF','par1','par2')

f=figure,
tiledlayout(2,1)
nexttile
plot(rangeLen,par1)
ylabel('\nu')
nexttile
plot(rangeLen,par2)
ylabel('\lambda')
xlabel('bar length')
saveas(f,'pcc3.png')
save('pcc3.mat','rangeLen','par1','par2')
%% MP set up:
foldSynth = 'barcolinull';
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
mkdir(strcat('output',timestamp));
if ispc % for PC get from the initialization
    sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP' ;
    SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
else %
%     sets.SCAMP_LINE = '/home/albyback/postdocData/test_transloc/SCAMP/';
    sets.SCAMP_LINE = '~/SCAMP/';

    SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
end

MIN_OVERLAP_PIXELS = minOverlap;
%%
k=1;
par1=[];
par2=[];

rLen = 1200;
nmpx=110;
sF = {0.95:0.01:1.05};
% for nmpx=rangePSF % 3 loops each for parameter
% for sF=rangeSF % 3 loops each for parameter
for rLen = rangeLen
    PSF_WIDTH_PIXELS = 300/nmpx;
    RAND_LENGTH_MIN= rLen;
    [~,barcodeGen] = gen_random(NUM_RAND_FRAGMENTS,PSF_WIDTH_PIXELS,RAND_LENGTH_MIN,0);
    %
    barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);
    % have to save separately..
    [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF{1},foldSynth);
    % skip one barcode
    [names2, baridx2] = Core.save_long_bar_txt_skip(barStruct,1:length(barStruct),foldSynth);


    % delete(gcp('nocreate'))
    numWorkers = 32;
    tic
    import Core.calc_overlap;
    [mpI1,mp1,maxMP] = calc_overlap(barStruct,timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,names2,namesBar);
    toc
    
    
    for ii=1:length(baridx2)
    baridx2{ii} = baridx2{ii}(1:end-MIN_OVERLAP_PIXELS+1);
    end
    % we want to create a nicer structure for all-to-all comparison and
    % contains easily accesible data.
%     tic
    import Core.mp_res_to_struct;
    overlapStruct = mp_res_to_struct(mp1,mpI1,baridx2,stridx,MIN_OVERLAP_PIXELS,sF{1},barStruct);
%     toc
     
    PCC_OVERLAP = reshape([overlapStruct.score], size(overlapStruct,1),size(overlapStruct,2));
    overlaplen = reshape([overlapStruct.overlaplen], size(overlapStruct,1),size(overlapStruct,2));
    bestBarStretch = reshape([overlapStruct.bestBarStretch], size(overlapStruct,1),size(overlapStruct,2));
    
    % lenOverlap = 
    
    import Zeromodel.beta_ev_params;
    
    % fit the parameters:
    maxpcc = PCC_OVERLAP(~isnan(PCC_OVERLAP));
    [parameters] = beta_ev_params(maxpcc, minOverlap/3);

    [parameters(1)/minOverlap  parameters(2)]
    par1(k) = parameters(1)/minOverlap;
    par2(k) = parameters(2);
    k=k+1;
end
%%
f=figure,
tiledlayout(2,1)
nexttile
plot(rangePSF,par1)
ylabel('\nu')
nexttile
plot(rangePSF,par2)
ylabel('\lambda')
xlabel('nm/px')
saveas(f,'mp1.png')
save('mp1.mat','rangePSF','par1','par2')


f=figure,
tiledlayout(2,1)
nexttile
plot(par1)
ylabel('\nu')
nexttile
plot(par2)
ylabel('\lambda')
xlabel('sF')
saveas(f,'mp2.png')
save('mp2.mat','rangeSF','par1','par2')


f=figure,
tiledlayout(2,1)
nexttile
plot(rangeLen,par1)
ylabel('\nu')
nexttile
plot(rangeLen,par2)
ylabel('\lambda')
xlabel('bar length')
saveas(f,'mp3.png')
save('mp3.mat','rangeLen','par1','par2')

%% plot
    [parameters] = beta_ev_params(maxpcc, minOverlap/3);

xx =0:0.01:1;
import Zeromodel.beta_ev_pdf
values = beta_ev_pdf(xx, parameters(1), 1, parameters(2));
figure,plot(xx,values)

hold on
histogram(maxpcc,'Normalization','pdf')