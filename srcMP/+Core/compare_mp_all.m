function [mpI1individual, mp1individual, maxMP, stridxindividual, compStr,calcLengths] = compare_mp_all(theoryStruct,barcodeGen,minLen,ix,timestamp,sF,MIN_OVERLAP_PIXELS,numWorkers)
    % compare_mp_all
    %   Compare all barcodes against the theory.
    %   If barcode length is less than minLen, don't calcularw

    %   Args:
    %       theoryStruct,barcodeGen,minLen,ix,timestamp,sF,MIN_OVERLAP_PIXELS,numWorkers
    %   Returns:
    %
    %       mpI1individual,mp1individual,maxMP,stridxindividual,compStr

    barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    calcLengths = find(barLens>minLen);
    barcodeGen = barcodeGen(barLens>minLen);
    lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
    
    barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
        cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

    SCAMP_LINE = 'SCAMP';
    
    
    foldSynth = 'barcoli';

    % [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);
    [namesBar, stridx, baridx] = Core.save_bars_rescaled_txt_stack(barStruct,sF,foldSynth);
    stridx = stridx(1:end-MIN_OVERLAP_PIXELS+1);
    baridx = baridx(1:end-MIN_OVERLAP_PIXELS+1);
    
    % names2{1} = ;// could also stack the two theories (esp only if we're
    % interested in best match)
    names2 = arrayfun(@(x) theoryStruct{ix}.filename, 1:length(namesBar),'un',false);
    % names2 = arrayfun(@(x) theoryStruct{x}.filename, 1:length(theoryStruct),'un',false);
    
    tic
    import Core.calc_overlap;
    [mpI1,mp1,maxMP] = calc_overlap([],timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,namesBar,names2);
    toc

    sFList = [1:length(sF) -(1:length(sF))];
    
    % % %% Check the re-scaling factors along the barcode. 
    strFac = cell(1,length(barStruct));
    pos = cell(1,length(barStruct));
    bestSTF = zeros(1,length(barStruct));
    compStr = cell(1,length(barStruct));
    mp1individual = cell(1,length(barStruct));
    mpI1individual = cell(1,length(barStruct));
    stridxindividual = cell(1,length(barStruct));
    for i=1:length(barStruct)
        mp1individual{i} = mp1{1}(baridx==i);
        mpI1individual{i} = mpI1{1}(baridx==i);
        stridxindividual{i} = stridx(baridx==i);
    
        locs = [];
        meanMP=[];
        maxMp = [];
        for j=1:length(sFList) % length re-scaling based on average of a few highest values.
            mpS = mp1individual{i}(stridxindividual{i}==sFList(j));
            [mpS,locs{j}] = sort(mpS,'descend','MissingPlacement','last');
            meanMP(j) = nanmean(mpS(1:min(end,20))); %1 would be standard
            maxMp(j) = mpS(1);
            locs{j} = locs{j}(~isnan(mpS));
        end
        compStr{i}.meanMP = meanMP;
        [aval,bstr] = max(compStr{i}.meanMP); % find best re-scale factor
    
        strFac{i} = sFList(bstr); % +1??
    
        pos{i} = locs{bstr}(1);
    
        mpIs = mpI1individual{i}(stridxindividual{i}==sFList(bstr));
    
        bestSTF(i)=  sF(abs(strFac{i}(1)));
        compStr{i}.pB = mpIs(locs{bstr}(1))'+1;
    
    
        compStr{i}.pA =  pos{i};% position on strFac - from this can extract overlap pos 
        pos1 = find(barStruct(i).rawBitmask==1,1,'first');

        compStr{i}.pos = mpIs(locs{bstr}(1))'+1-compStr{i}.pA+1-pos1; % position w.r.t. unscaled barcode
    
        compStr{i}.or = sign( strFac{i}); % should be able to plot exact match..
        compStr{i}.idx = ix;
        compStr{i}.bestBarStretch =   bestSTF(i);
        compStr{i}.lengthMatch =  lengths(i)*bestSTF(i);
        compStr{i}.maxcoef = maxMp(bstr)';
        compStr{i}.w = MIN_OVERLAP_PIXELS;
        compStr{i}.allposA = locs{bstr};
        compStr{i}.allposB = mpIs(locs{bstr})+1;
    end

% % %% Check the re-scaling factors along the barcode. 
% % strFac = cell(1,length(mp1));
% % pos = cell(1,length(mp1));
% % bestSTF = zeros(1,length(mp1));
% % 
% % for i=1:length(mp1)
% %     [a,b] = sort(mp1{i},'Desc');
% %     strFac{i} = stridx{i}(b(1:(sum(stridx{i}==1)-MIN_OVERLAP_PIXELS )));
% %     pos{i} = mpI1{i}(b(1:(sum(stridx{i}==1)-MIN_OVERLAP_PIXELS )));
% %     bestSTF(i)= sF(abs(strFac{i}(1)));
% % end
% %     
% % % check if all position along the genome have the same re-scaling factor. 
% % % MP can be extended to full overlap, and we find best re-scaling factor
% % % for all those positions
% % % The same procedure can be done for the all to all comparison.    
% % 
% % 
% % %% for experiment
% % [xId,yId] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
% % %%
% % % we compare xId with yId. yId is un-rescaled and xId is re-scaled
% % % xId = 195
% % % yId = 163
% % % we find the positions on mp1 and mpI1
% % positionsBar = find(baridx2{xId}==yId);
% % 
% % mpBar = mp1{xId}(positionsBar);
% % mpIBar = mpI1{xId}(positionsBar);
% % rescaleBar = mpI1{xId}(positionsBar);
% % 
% % mpIBar = mpIBar(mpIBar>-1);
% % 
% % 
% % % figure2:
% % strFactors = stridx{xId}(mpIBar+1);
% % 
% % 
% % % figure 1: mpI positions
% % % figure,plot(mpIBar)
% % % 
% % % import Core.get_two_bars_alignment_params;
% % % [bar1Name, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{xId}, mpI1{xId}, ps, baridx2{xId},sF);
% % 
% % % should only keep the re-scaling factors for the pixels where the barcodes
% % % are overlapping
% %  sortedidx2 = find(pksUnique1{xId}==yId);
% % %     sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);
% % 
% % pos = LOCS1{xId}(pksUniquePos1{xId}(sortedidx2)); % position on MP.
% % import Core.get_two_bars_alignment_params;
% % 
% % [~, pA, pB,rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{xId},mpI1{xId},pos, baridx2{xId},sF); 
% % % pos1 = -pA+1:-pA+length(bRescaled);
% % % pos2 = -pB+1:-pB+length(b2);
% % 
% % 
% % %
% % 
% % vec = [yId xId];
% % % [bar1idx,bar2idx] = [29],[49];
% % bar1idx=vec(2);
% % bar2idx=vec(1);
% % sortedidx2 = find(pksUnique1{bar1idx}==bar2idx);
% % import Core.plot_best_match;
% % [f,pos1, bar1, pos2, bar2,pcc1,pcc2] = plot_best_match(sortedidx2,stridx,mpI1{bar1idx},LOCS1{bar1idx},pksUniquePos1{bar1idx},bar1idx,baridx2{bar1idx},sF,barStruct,MIN_OVERLAP_PIXELS);
% % pcc1
% % pcc2
% % 
% % [C,IA,IB] = intersect(pos1,pos2);
% % xlim([C(1) C(end)])
% % 
% % nexttile([1 2]),plot(-pB:-pB-1+length(strFactors),[strFactors])
% % xlim([C(1) C(end)])
% % nexttile([1 2]),plot(-pB:-pB-1+length(strFactors),[mpBar(~isnan(mpBar))])
% % xlim([C(1) C(end)])
% % 
% % % xlim([-pA+1 -pA+length(bar1)])
% % 
% % % todo: estimate PSF/NMBP/NMPX based on data
% % 

end

