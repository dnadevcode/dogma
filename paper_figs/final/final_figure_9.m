% validation of fig_6


% First generate theory:
w = 300;% minimum overlap length
sF = 0.9:0.025:1.1;
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

fastas = {'018_final_polish.fasta','DA32087.fasta','DA68335.fasta'};
theoryIdxs = {1, 3, 2, nan, 2 , 2, 3}; % known theory indexes
dataSetIdx = 2; % 5,... test datasets (part of the data)

nmPerPx = 110;
nmbp = 0.2; % ? 

% quickly calculate theory for given single fasta file
[theoryStructRev,theoryStruct,bT] = prep_thry_for_local_comp(fastas(theoryIdxs{dataSetIdx}), nmbp, nmPerPx, 1);


import Core.update_sf_barset;

%% Validation. 
theoryStruct.rawBitmask = bT{1}.rawBitmask;
% tS.length = length(bT{1}.rawBitmask);
pDif = cell(1,length(barIslands));
sorti =  cell(1,length(barIslands));
val =  cell(1,length(barIslands));
for idx = 1:length(barIslands)

    %%1) Run comparison (all barcodes)
    bars = barcodeGen(barIslands{idx}); 
    barx = barIslands{idx};
    sets.comparisonMethod = 'mass_pcc';
    sets.nuF = 0.08;
    [compI,rezI,~] = compare_to_t(bars,theoryStruct,sF,sets); % maybe calculate Stouffer score here also?


%     N = 1000;
% 
% import Validation.calc_distance_scaled; % need to properly scale things when calculating distances for validation
% posShift = calc_distance_scaled(oS,sortedIdsGood(1:N), compI,theoryStruct.length, N);
%     
%     compI = compI(barIslands{idx});


    [a] = cellfun(@(x) x.pval,compI);
    [minV,minPos] = min(a);
    %% 
    selBarId = minPos;

    %% extract tables
    bbS = compI{selBarId}.bestBarStretch;  % best bar stretch (bbS)
    bbO = (compI{selBarId}.or(1)~=barcodeIslandsData{idx}(selBarId,3))+1;

    [bbSAll] = update_sf_barset(barcodeIslandsData{idx}, bbS/barcodeIslandsData{idx}(selBarId,4), bbO);

    % need to also update positions for all the synCur barcodes based on
    import Core.synth_to_table;
    [tableS] = synth_to_table(compI);

    % update synth table. Want to have the same
    [tableSUpd] = update_sf_barset(tableS, tableS(selBarId,4)/bbS, 1);


    bbSAll(:,1:2) = bbSAll(:,1:2)-bbSAll(selBarId,1)+tableSUpd(selBarId,1); %
   
    bestPosFound = bbSAll(:,1);
    
    posTrue = tableSUpd(:,1);

    pDif{idx} =  min([sqrt((posTrue-bestPosFound).^2) sqrt((posTrue+theoryStruct.length-bestPosFound).^2)...
        sqrt((posTrue-theoryStruct.length-bestPosFound).^2)]');

    [val{idx},sorti{idx}] = sort(posTrue);

end


%%

% scaled table to map:
import Core.consensus_from_table
[outConsensus] = consensus_from_table(bbSAll,bars);

[pos,idxv] = sort(arrayfun(@(x) find(~isnan(outConsensus(x,:)),1,'first') +(find(~isnan(outConsensus(x,:)),1,'last')-find(~isnan(outConsensus(x,:)),1,'first'))/2,1:size(outConsensus,1)));
% 
consensusToPlot1 = outConsensus(idxv,:);

imagesc(consensusToPlot1);colormap(gray)


foundPos = bestPosFound(selBarId);

shift = -min(bbSAll(:,1))+1;

subMat = consensusToPlot1(:,shift+foundPos:shift+foundPos+bbSAll(selBarId,2)-bbSAll(selBarId,1));

% B=subMat(sum(isnan(subMat),2));

% any(isnan(subMat'))

subMat(all(isnan(subMat), 2), :) = [];
% subMat(any(isnan(subMat), 2), :) = [];
%%
f = figure;
tiledlayout(5,1)
nexttile
imagesc(tableS(selBarId,1):tableS(selBarId,2),1,theoryStruct.rawBarcode(tableS(selBarId,1):tableS(selBarId,2)));colormap(gray)
axis off
title('(A) Part of theory corresponding to optimally matching barcode','Interpreter','latex')
nexttile
subMat(subMat<-4) = nan;
imagesc(tableS(selBarId,1):tableS(selBarId,2),1,subMat);colormap gray
title('(B) Corresponding part of block representation of barcode island','Interpreter','latex')

nm = {'(C)','(D)','(E)','(F)'}

for idx=1:length(val)
    nexttile
    plot(val{idx},     pDif{idx}(sorti{idx}),'blackx')
    xlabel('Position (px)')
    ylabel('pd')

    title([nm{idx} ' Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')

end
% if bbSAll
    print('FIGS/Fig9.eps','-depsc','-r300');


%%
figure,plot(posTrue,bbSAll(:,1),'x');xlabel('ground truth (px)');ylabel('placement position (px)');title('Dot plot for barcode island (2)')

% Local similarity
wSim = 500;
[MPthr, MPIthr1] = local_compare_thry(wSim,[], nmbp, nmPerPx, 1,32,theoryStructRev)

mpVals = MPthr{1}(1:theoryStruct.length-wSim+1);
% mpVals = MPthr{1}(1:theoryStruct.length);


figure
% tiledlayout(2,1)
% nexttile
% idx =2
% plot(val{idx},     pDif{idx}(sorti{idx}),'blackx')
% xlabel('Position (px)')
% ylabel('pd')
% 
% title([' Distance for barcode island (',num2str(idx),')'], 'Interpreter','latex')
% xlim([1 length(theoryStruct.rawBarcode)])

% nexttile
hold on
title('Local similarity along the theory barcode')
plot(mpVals) % depends if we flip theory to get 
xlabel('Position (px)')

xlim([1 length(theoryStruct.rawBarcode)])

