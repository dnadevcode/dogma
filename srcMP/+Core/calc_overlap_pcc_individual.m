function [oS,overlapStruct] = calc_overlap_pcc_individual(barcodeGen, barcodeGen2, sF, minOverlap)
    % extracting scores for different length overlaps for pcc

    % version that sorts based on individual p-values.
    
import SignalRegistration.masked_pcc_corr;
% barcodeGen2 = barcodeGen;

% minOverlap = 300;
% overlapStruct = [];%cell(1,length(mp1));
overlapStruct = cell(length(barcodeGen),length(barcodeGen2));%cell(1,length(mp1));

psfPar = 0.03; % depends on psf. so possibly adjust for the re-scale factor too
f = @(x,sigma) 1/2.*(1+erf(x./(sqrt(2).*sigma)));

for i=1:length(barcodeGen)
    i
    A = arrayfun(@(y) imresize(barcodeGen{i}.rawBarcode,'Scale' ,[1 y]),sF,'un',false);
    B = arrayfun(@(y) imresize(barcodeGen{i}.rawBitmask,'Scale' ,[1 y]),sF,'un',false);
%     i
    for j=1:length(barcodeGen2)
%         j

%         if i~=j
%             scores = nan(1,ceil(max(cellfun(@(x) sum(x.rawBitmask),barcodeGen))*sF(end)));
%             nE = zeros(1,ceil(max(cellfun(@(x) sum(x.rawBitmask),barcodeGen))*sF(end)));

            barB = barcodeGen2{j}.rawBarcode;
            wB = barcodeGen2{j}.rawBitmask;
            scoreCur = zeros(1,length(A));
            pccCur = zeros(1,length(A));
            orCur = zeros(1,length(A));
            posCur =  zeros(1,length(A));
            numElts =cell(1,length(A));
            pos = zeros(1,length(A));
            parfor s = 1:length(A)
                extraL = length(A{s});%max(length(A{s})-minOverlap,length(A{s})-length(barB));
                barB2 = [zeros(1,extraL) barB];
                wB2 = [zeros(1,extraL) wB];
                
%                 tic
                [ xcorrs, numElts{s} ] = masked_pcc_corr( A{s},barB2,B{s},wB2,minOverlap ); %todo: include division to k to reduce mem
%                 toc
%                 xcorrs(numElts{s}>4)=nan;
%                 tic % es
%                 f = @(x,n) (1-x.^2).^((n-4)/2)./beta(1/2,1/2*(n-2));
%                 scaledXcorrs = xcorrs./log(numElts{s});
%                 f = @(x,n) 1-1/2*(1+erf(x/(sqrt(2)*n)));
                scaledXcorrs = f(xcorrs,1./sqrt(psfPar*numElts{s}-2) );

%                 toc
                
                % number of elements and maximum for each overlap length
%                 tic
%                 allelts = xcorrs(:);
%                 [valE, posE]= sort(numElts{s}(:));
%                 [c1,ia1,ic1] = unique(valE);
%                 numECur = diff(ia1);
%                 
%                 nonNans = c1(~isnan(c1));
%                 nE(nonNans) = nE(nonNans)+numECur(1:length(nonNans))';
%                 mc1 =min(c1);
%                 maxc1 = max(c1);
%                 for idd=mc1:maxc1
%                     scores(c1(idd-mc1+1)) = nanmax([scores(c1(idd-mc1+1)) allelts(posE(ia1(idd-mc1+1):ia1(idd-mc1+2)-1))']);
%                 end
%                 toc
        
%                 [lenOverlap numE] = unique(numElts{s});
%                 tic
                % here xcorrs should be sorted.
                [a,b] = max(scaledXcorrs);
%                 tic
                [c,d] = max(a); 
                orCur(s) = b(d);
                posCur(s) = d-extraL;
                pos(s) = d;
                scoreCur(s) = c;%xcorrs(b(d),d); % where to save c?
                pccCur(s) = xcorrs(b(d),d);
%                 toc
    
            end
            [maxS,posMax] = max(scoreCur);
            overlapStruct{i,j}.score = maxS;
            overlapStruct{i,j}.fullscore = maxS;

            overlapStruct{i,j}.pB = posCur(posMax); % position on root (un-rescaled barcode)
            overlapStruct{i,j}.pos = overlapStruct{i,j}.pB;
            overlapStruct{i,j}.idx = j;
            overlapStruct{i,j}.pA = 1; % position on root (un-rescaled barcode)

            overlapStruct{i,j}.or = orCur(posMax); %pos on rescaled bar
            overlapStruct{i,j}.bestBarStretch = sF(posMax); %orientation
            overlapStruct{i,j}.overlaplen = numElts{posMax}(orCur(posMax),pos(posMax));
            overlapStruct{i,j}.lengthMatch =    overlapStruct{i,j}.overlaplen;
            overlapStruct{i,j}.h = numElts{posMax}(orCur(posMax),pos(posMax));
            overlapStruct{i,j}.lenA  = sum(B{posMax});
            overlapStruct{i,j}.lenB  = sum(wB);
%             overlapStruct{i,j}.xcorrs = xcorrs;
%             overlapStruct{i,j}.numElts = numElts;
            overlapStruct{i,j}.sc = pccCur(posMax);
            
%             overlapStruct{i,j}.indScores = scores;
%             overlapStruct{i,j}.indOverlaps = nE;


%         else
%             overlapStruct{i,j}.score  = nan;
%             overlapStruct{i,j}.fullscore  = nan;
%             overlapStruct{i,j}.overlaplen  = nan;
%             overlapStruct{i,j}.pA  = nan;
%             overlapStruct{i,j}.or  = nan;
%             overlapStruct{i,j}.pB  = nan;
% %             overlapStruct{i,j}.score  = nan;
%             overlapStruct{i,j}.lenA  = nan;
%             overlapStruct{i,j}.lenB  = nan;
%             overlapStruct{i,j}.bestBarStretch = nan;
%               overlapStruct{i,j}.xcorrs = nan;
%             overlapStruct{i,j}.numElts = nan;
%               overlapStruct{i,j}.sc = nan;
%             overlapStruct{i,j}.indScores = nan(1,ceil(max(cellfun(@(x) sum(x.rawBitmask),barcodeGen))*sF(end)));
%             overlapStruct{i,j}.indOverlaps = nan(1,ceil(max(cellfun(@(x) sum(x.rawBitmask),barcodeGen))*sF(end)));
%         end

    end
end

oS = [];
for i=1:size(overlapStruct,1)
    for j=1:size(overlapStruct,2)
        oS(i,j).score  = overlapStruct{i,j}.score;
        oS(i,j).fullscore =  overlapStruct{i,j}.fullscore;
        oS(i,j).overlaplen  = overlapStruct{i,j}.overlaplen;
        oS(i,j).pA  = overlapStruct{i,j}.pA;
        oS(i,j).or  = overlapStruct{i,j}.or;
        oS(i,j).pB  = overlapStruct{i,j}.pB;
        %             oS(i,j).score  = overlapStruct{i,j}.score;
        oS(i,j).lenA  = overlapStruct{i,j}.lenA;
        oS(i,j).lenB  = overlapStruct{i,j}.lenB;
        oS(i,j).bestBarStretch  = overlapStruct{i,j}.bestBarStretch;
%                oS(i,j).xcorrs=   overlapStruct{i,j}.xcorrs ;
%             oS(i,j).numElts =  overlapStruct{i,j}.numElts;
           oS(i,j).sc = overlapStruct{i,j}.sc;
%         oS(i,j).indScores  = overlapStruct{i,j}.indScores;
%         oS(i,j).indOverlaps  = overlapStruct{i,j}.indOverlaps;

    end
end

% end
% 
% 
% 
%            subMP = mp1{k}(baridx2{k}==iy);
%             subMPI = mpI1{k}(baridx2{k}==iy);
%             [mpS,locs] = sort(subMP,'descend','MissingPlacement','last');
%             overlapStruct(k,iy).score = mpS(1);
%             overlapStruct(k,iy).pB = locs(1); % position on root (un-rescaled barcode)
% 
%             %
%             srtIdx =  stridx{k}(subMPI(locs(1))+1);
%             posfirst = find(stridx{k}==srtIdx,1,'first'); % if we take just the best
% 
%             overlapStruct(k,iy).pA = subMPI(locs(1))+1-posfirst+1;
%             overlapStruct(k,iy).or = sign(srtIdx);
%             overlapStruct(k,iy).bestBarStretch = sF(abs(srtIdx));
%             overlapStruct(k,iy).overlaplen = MIN_OVERLAP_PIXELS;
% 
% 
%         else
%             overlapStruct(k,iy).score  = nan;
%             overlapStruct(k,iy).fullscore  = nan;
%             overlapStruct(k,iy).overlaplen  = nan;
%             overlapStruct(k,iy).pA  = nan;
%             overlapStruct(k,iy).or  = nan;
%             overlapStruct(k,iy).pB  = nan;
%             overlapStruct(k,iy).score  = nan;
%             overlapStruct(k,iy).lenA  = nan;
%             overlapStruct(k,iy).lenB  = nan;
% 
%         end
%     end
% end
% 
