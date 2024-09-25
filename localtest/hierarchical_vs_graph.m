
% first generate some random data
NUM_RAND_FRAGMENTS = 100; % number of random fragments
NUM_RAND_FRAGMENTS*(1-NUM_RAND_FRAGMENTS/10000)^500

[barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr] = gen_rand(NUM_RAND_FRAGMENTS);


barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen','un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen','un',false)]',{'rawBarcode','rawBitmask'},2);

%% visualization/testing

bpPx = 500;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], synthStr, [], [], [],theoryStruct{1}.length,1);%1/bpPx*10^6
% 
psc = reshape([synthStr2(:).score], size(synthStr2,1),size(synthStr2,2)); % from bg_test_1
filtM = triu(psc);
filtM(filtM==0) = nan;
pscores = filtM(:);

% 
% create barset island
import Core.create_barset;
[barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = create_barset(pscores,sum(~isnan(filtM(:))), barcodeGen, synthStr2,synthStr);

nonEmpty= find(cellfun(@(x) ~isempty(x),barcodeIslands));

% plot of detected islands
import Plot.islandsPosStruct;
posStruct = islandsPosStruct(barcodeIslandsData(nonEmpty),barcodeIslands(nonEmpty));
fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(10000,1,length(barcodeIslandsData(nonEmpty)))),1);%1/bpPx*10^6

% general positions of these
import Plot.plot_best_pos;
fig1 = plot_best_pos([], synthStr(cell2mat(barcodeIslands)), [], [], [],10000,1);%1/bpPx*10^6

% import Plot.
%      G=create_graph(barcodeGen,PCC_OVERLAP,sortedIds,NN);

%% TODO: align barcode positions to compare detected vs true pos
badIdx = find(cellfun(@(x) ~isempty(x),badDataInfo));
ix = 12;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], synthStr(cell2mat(barIslands{ix}.barcodeIslands)), [], [], [],10000,1);%1/bpPx*10^6

% posStruct = islandsPosStruct(barcodeIslandsData,barcodeIslands);
fig1 = plot_best_pos([], barIslands{ix}.posStruct, [], [], [],cumsum(repmat(2000,1,length(barcodeIslandsData))),1);%1/bpPx*10^6


% locs = find(cellfun(@(x) ~isempty(x),barcodeIslandsData));
% figure;hold on
% for i=1:size(barcodeIslandsData{locs(end)},1)
%     plot([barcodeIslandsData{locs(end)}(i,1) barcodeIslandsData{locs(end)}(i,2)],[i i],'red-')
% end

% %
% i = 13;
% j = 16;
import Core.plot_match_pcc;
% [f] = plot_match_pcc(barStruct, synthStr2,i,j,barStruct);
% %%

% % plot
idd =[6 7];
i = idd(1);
j = idd(2);
% [outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen(idd),synthStr(idd),theoryStruct,datestr(clock(), 'yyyy-mm-dd_HH_MM_SS'));
plot_match_pcc(barStruct, synthStr2,idd(1),idd(2),barStruct);
% figure,plot(outConsensus')
%%
% zscore HC
% % tic
res = compute_bar_pairs_zscore(barcodeGen);
% % toc

% PCC/MP graph
[res2] = compute_bar_pairs_graph(barcodeGen);
% 
function [barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr] = gen_rand(NUM_RAND_FRAGMENTS);

% 10/01/22 Synthetic
% version of chrom_assembly_script_vs_ref_synthetic
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% STEP1: gen synthetic data

%% gen synthetic
import Nullmodel.gen_random_fragments;

% parameters
TOTAL_RAND_LENGTH = 10000; % total length
PSF_WIDTH_PIXELS = 4.7; % psf
% NUM_RAND_FRAGMENTS = 20; % number of random fragments
MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
FRAGMENT_LENGTH_STD = 100; %std of length of fragment
ADDED_NOISE_MEAN = 0.4; % additive noise mean
ADDED_NOISE_STD = 0.1; % additive noise std
FRAGMENT_STRETCH_STD = 0;%0.02; % fragment length-rescale factor std 0.02
IS_CIRCULAR = 1; % whether circular barcode
    
% generate barcodes /also add pcc?
[barcodeGen, refBarcode, origPos, origFlip,origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD,IS_CIRCULAR);


% plot_quick(6,refBarcode,barcodeGen,origPos, origFlip,1./origStr)

% only keep barcodes longer than minimum length.
% minLen = 300; %150kb? or less? depends on application
% barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
% barcodeGen = barcodeGen(barLens>minLen);
% lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);


% plot the fragment positions. Requires hca function
synthStr = cell(1,length(barcodeGen));
for i=1:length(barcodeGen)
    synthStr{i}.idx = 1;
    synthStr{i}.pos = origPos(i);
    synthStr{i}.maxcoef = 0.99; % calculate max coef against thry
    synthStr{i}.lengthMatch = round(length(barcodeGen{i}.rawBarcode)/origStr(i));
    synthStr{i}.lengthUnrescaled = length(barcodeGen{i}.rawBarcode);
    synthStr{i}.bestBarStretch = 1./origStr(i);
    synthStr{i}.or = double(origFlip(i))+1;
    synthStr{i}.rf = origStr(i);
end

% if circ, check if loop around
synthStrLin = synthStr;
% for i=1:length(barcodeGen)
%     if synthStr{i}.pos+synthStr{i}.lengthMatch-1>TOTAL_RAND_LENGTH
%         synthStrLin{i}.pos = synthStrLin{i}.pos-TOTAL_RAND_LENGTH;
%     end
% end
% expected overlap positions
%
synthStr2 =[];
for i=1:length(barcodeGen)
    for j=1:length(barcodeGen)
        if i~=j
            synStemp = synthStrLin;
            % make sure that synStemp{i}.pos-synStemp{j}.pos close to each
            % other
            
            % if distance greater than total length, move the barcode which
            % is more to the right, to the left
            if synStemp{i}.pos>synStemp{j}.pos
                idx = i;
            else
                idx = j;
            end
            if abs(synStemp{i}.pos-synStemp{j}.pos)>abs(abs(synStemp{i}.pos-synStemp{j}.pos)-TOTAL_RAND_LENGTH)
                synStemp{idx}.pos = synStemp{idx}.pos-TOTAL_RAND_LENGTH;
            end

                 
            rF1 = (synthStr{j}.rf); % rescale A to B length
            rF2 = 1;%1/origStr(j);
            
% 

            % length with proper length re-scaling where B is not
            % re-scaled and A re-scaled to B
            lpA =  round(synStemp{i}.lengthMatch*rF1);
            lpB =  round(synthStr{j}.lengthUnrescaled);
         
            synthStr2(i,j).posB = 1;
            synthStr2(i,j).posA = (synStemp{i}.pos-synStemp{j}.pos);

            if synthStr{j}.or==2 % should have j always with or=0 and then switch, or allow it to be flipped?
                pT  = (synStemp{j}.pos+synStemp{j}.lengthMatch-synStemp{i}.pos-synStemp{i}.lengthMatch)+1; % circular?
                synStemp{i}.or =  -synStemp{i}.or+3;
            else
                pT = synthStr2(i,j).posA+1; % 0 should mean identical position
            end
%                 synthStr2(i,j).posA = round((pT)*origStr(j)); % start position is adjusted by origStr(j)
%             else
%%
                % depends if pos or neg.
            if pT > 0
                synthStr2(i,j).posA =round(pT*origStr(j)); % start position is adjusted by origStr(j) to get correct scaling
            else
                synthStr2(i,j).posA =round(pT*origStr(j)); % correct by 1 if before was from 0 
            end
            
          % consider the jth barcode as always at 1 (posB)

            
            pA =  synthStr2(i,j).posA ;
            pB =  synthStr2(i,j).posB ;
            
            st = max(pA,pB)-pB+1;%-( pA-pB+1); % start along un-rescaled reference
            stop = min(pA+lpA-1,pB+lpB-1)-pB+1;
            
      
            re = barcodeGen{i}.rawBarcode;
            re = interp1(re, linspace(1, length(re), round(length(re)*rF1*synthStr{i}.bestBarStretch)));
%             reb = barcodeGen{i}.rawBitmask;
%             reb = interp1(reb, linspace(1, length(reb), round(length(reb)*rF1)));
            rt = barcodeGen{j}.rawBarcode;
            rt = interp1(rt, linspace(1, length(rt), round(length(rt)*rF2)));

            a = [   re re];
            b =[ rt rt];
            if synStemp{i}.or==2
                a = fliplr(a);
            end
%             if synthStr{j}.or==1 % should have j always with or=0 and then switch, or allow it to be flipped?
%                 b = fliplr(b);
%             end
            aFul = a(st-synthStr2(i,j).posA+1:stop-synthStr2(i,j).posA+1);
            bFul = b(st:stop);
            
            % position of overlap
            synthStr2(i,j).pA = st-synthStr2(i,j).posA+1 ;
            synthStr2(i,j).pB = st ;
% 
%             figure,plot(zscore(aFul))
%             hold on
%             plot(zscore(bFul))
%             legend({'i','j'})
            %%
            
%             figure,plot(-synthStr2(i,j).pA:-synthStr2(i,j).pA+length(re)-1,re)
%             hold on
%             plot(-synthStr2(i,j).pB:-synthStr2(i,j).pB+length(rt)-1,rt)

% %       
            synthStr2(i,j).fulloverlapPosA = pB-st+1:pB+stop-1;
            synthStr2(i,j).fulloverlapPosARescaled = pA-st+1:pA+stop-1;
            synthStr2(i,j).fullscore = zscore(aFul(10:end-9),1)*zscore(bFul(10:end-9),1)'/length(aFul(10:end-9));

            synthStr2(i,j).overlaplen = length(aFul);
            if synthStr2(i,j).overlaplen > 150
                synthStr2(i,j).score =   synthStr2(i,j).fullscore ;
            else
               synthStr2(i,j).score = nan;
            end
            synthStr2(i,j).bestBarStretch = rF1*synthStr{i}.bestBarStretch;
            synthStr2(i,j).or = synStemp{i}.or;
            
            synthStr2(i,j).lenB = length(barcodeGen{j}.rawBarcode);
            synthStr2(i,j).lenA = length(interp1(barcodeGen{i}.rawBarcode, linspace(1, length(barcodeGen{i}.rawBarcode), round(length(barcodeGen{i}.rawBarcode)*synthStr2(i,j).bestBarStretch))));
%             synthStr2(i,j).fullscore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
%             overlapStruct(k,iy).overlaplen = length(aFul);
%             overlapPos = 
        else
            synthStr2(i,j).fullscore =nan;
            synthStr2(i,j).score =nan;

        end   
    end
end
% nd

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

end


function plot_quick(ix,refBarcode,barcodeGen,origPos, origFlip,origStr)
%     ix =1;
    re = barcodeGen{ix}.rawBarcode;
    re = interp1(re, linspace(1, length(re), round(length(re)*origStr(ix))));
    refBarcode = [refBarcode refBarcode];
    rt = refBarcode(origPos(ix):origPos(ix)+length(re)-1);
    if origFlip(ix)==1
        re = fliplr(re);    
    end
    
    figure,plot(rt)
    hold on
    plot(re)
    
    zscore(re(10:end-9))*zscore(rt(10:end-9)')/length(rt(10:end-9))
end


function [res] = compute_bar_pairs_zscore(barcodeGen)
    %   compute_bar_pairs_zscore

    %%
    MIN_OVERLAP_PIXELS = 150; % minimum number of overlap between two barcodes
    STRETCH_MAX = 0; % maximum stretch
    STRETCH_STEP = 0.01; % stretch step
    MIN_BARCODE_SCORE_THRESHOLD = 3;
%     MIN_BARGROUP_SCORE_THRESHOLD = 7; % threshold for belonging to bargroup 
%     MIN_REFGROUP_SCORE_THRESHOLD = 0; % threshold for belonging to refgroup (so 0 means everything is mapped somewhere)
    P_COMBINE_METHOD = 'stouffer';
    PRINT_TO_WINDOW = 1;


    barcodes = cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);
    bitmasks = cellfun(@(x) x.rawBitmask,barcodeGen,'un',false);

% Mask outliers
% bitmasks = cellfun(@(x) [false(1, numUncert) true(1, size(x, 2) - 2 * numUncert) false(1, numUncert)], barcodes, 'un', 0);
% bitmasks = cellfun(@(x, y) y & abs(x - mean(x(y))) <= 2.5 * std(x(y)), barcodes, bitmasks, 'un', 0);

%%
    import Core.Bargroup
    % refId = length(barcodes)+1;
    % refBargroup = Bargroup(refId, 1, 0, 1, length(refBarcode));

    %% Run all to all comparison for all re-scaling factors


    [scoreRegistry,scoreRegistryIndex] = Core.compute_all_barcode_pair_scores([barcodes], [bitmasks], MIN_BARCODE_SCORE_THRESHOLD, STRETCH_MAX, STRETCH_STEP, MIN_OVERLAP_PIXELS, P_COMBINE_METHOD, PRINT_TO_WINDOW);

    res.scoreRegistry = scoreRegistry;
    res.scoreRegistryIndex = scoreRegistryIndex;
    
    [outputTree,~, mergeScores] = Core.bargroup_hc(barcodes, scoreRegistry, scoreRegistryIndex, 0, STRETCH_MAX, STRETCH_STEP, P_COMBINE_METHOD,MIN_OVERLAP_PIXELS, PRINT_TO_WINDOW);

    res.outputTree = outputTree;
    res.mergeScores = mergeScores;
end


function [res] = compute_bar_pairs_graph(barcodeGen)
    sF = 0.95:0.01:1.05;
    minOverlap = 150;
    scorethresh = 0.7; % instead use Stouffer like for HC
    tic
    import Core.calc_overlap_pcc_sort_m;
    [overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap);
    toc
    
    res.overlapStruct = overlapStruct;

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    mkdir(strcat('output',timestamp));
    
%     tic
    % [overlapStructMP] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap);
    [overlapStructMP] = calc_overlap_mp(barcodeGen',sF, minOverlap-50,timestamp)
%     toc

%     res.overlapStructMP = overlapStructMP;
    
    
    import Zeromodel.pval_for_struct;
    pscores = pval_for_struct(overlapStruct,scorethresh,0.05);
    
    import Nullmodel.local_sim_test;
    [curStd, pairs, scrs] = local_sim_test(pscores, barcodeGen', overlapStruct,minOverlap-50,150,3000,timestamp,overlapStructMP);

    psc = reshape([pscores], size(overlapStruct,1),size(overlapStruct,2)); % from bg_test_1

    filtM = triu(psc);
    filtM(filtM==0) = nan;

    
    import Core.create_barset;
    [barcodeIslandsData, barcodeIslands, badData] = create_barset(pscores,sum(~isnan(filtM(:))),barcodeGen,overlapStruct);


    locs = find(cellfun(@(x) ~isempty(x),barcodeIslandsData));
    figure;hold on
    for i=1:size(barcodeIslandsData{locs(end)},1)
        plot([barcodeIslandsData{locs(end)}(i,1) barcodeIslandsData{locs(end)}(i,2)],[i i],'red-')
    end

    import Core.create_overlap_graph
    [finalgraph,Ggraphs] = create_overlap_graph(pscores,sum(~isnan(filtM(:))),barcodeGen,curStd,100,overlapStruct);

    import Core.create_graph_coverage_map;
    [data,prev] = create_graph_coverage_map(finalgraph,overlapStruct,1)


    res.pscores = pscores;
    res.finalgraph = finalgraph;

end
% 
% function y=limfig
% imglimit=10;
% if length(findobj('type','figure'))<imglimit
%     y=figure; 
% 
% else
%     'too many figures already open'
%     return
% end
% end
