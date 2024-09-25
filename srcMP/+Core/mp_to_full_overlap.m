function [PCC_OVERLAP,len1,len2,lenOverlap,PCC_MP] = mp_to_full_overlap(stridx,pksUnique1,pksUniquePos1,LOCS1,mpI1,baridx2, MIN_OVERLAP_PIXELS, barStruct, sF)

%   mp_to_full_overlap - converts from local alignment to full overlap
%   alignment (i.e. PCC over all overlaping pixels)

% Returns:
%   PCC_OVERLAP - total overlap
%   len1 - length bar1
%   len2 - length bar2
%   lenOverlap - length overlap
%   PCC_MP

% ix,stridx,mpI1{ix1},LOCS1{ix1},pksUniquePos1{ix1},ix1,baridx2{ix1},sF,barStruct,MIN_OVERLAP_PIXELS

%% example
% barcode IDX_1 matches to unique barcodes in (unsorted) list pksUnique1{IDX_1}
% we find the index of IDX_2 in this list, 
% sIDX_2.
%
% This gives the index of barcode IDX_2 in index of PCC scores for IDX_1, 
% pksUniquePos1{IDX_1}(sIDX_2)
%
% which can then be converted to the position along the barcode
% LOCS1{IDX_1}(pksUniquePos1{IDX_1}(sIDX_2))
%
% we went extract the positions of barcodes with respect of each other
% import Core.get_two_bars_alignment_params;
% [bar1Name, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{IDX_1},mpI1{IDX_1}, pos, baridx2{IDX_1},sF);


PCC_OVERLAP = zeros(length(mpI1),length(mpI1));
len1 = zeros(length(mpI1),length(mpI1));
len2 = zeros(length(mpI1),length(mpI1));
lenOverlap = zeros(length(mpI1),length(mpI1));
PCC_MP = zeros(length(mpI1),length(mpI1));
for IDX_1=1:length(mpI1)
IDX_1
    for IDX_2 =1:length(mpI1)
        if IDX_1 ~= IDX_2
                sIDX_2 = find(pksUnique1{IDX_1}==IDX_2);

                if ~isempty(sIDX_2)
                    pos = LOCS1{IDX_1}(pksUniquePos1{IDX_1}(sIDX_2)); % position on MP.
                    % ts2 = IDX_2

                    import Core.get_two_bars_alignment_params;
                    [bar1Name, pA, pB, rescaleFactor, orSign] = get_two_bars_alignment_params(stridx{IDX_1},mpI1{IDX_1}, pos, baridx2{IDX_1},sF);

                    b2 = barStruct(bar1Name).rawBarcode(barStruct(bar1Name).rawBitmask);
                    b = barStruct(IDX_1).rawBarcode(barStruct(IDX_1).rawBitmask);

                    bRescaled = imresize(b,'Scale' ,[1 rescaleFactor]);

                    if orSign==-1
                        bRescaled = fliplr(bRescaled);
                    end
                        pos1 = -pA+1:-pA+length(bRescaled);
                        pos2 = -pB+1:-pB+length(b2);
                        bar2 = b2;

                        fP1 = find(pos1==1,1,'first');
                        fP2 = find(pos2==1,1,'first');
                    % figure; hold on
                    %     plot(zscore(bRescaled(fP1:fP1+MIN_OVERLAP_PIXELS-1)),'black')
                    %     plot(zscore(bar2(fP2:fP2+MIN_OVERLAP_PIXELS-1)),'red')

                        PCC_MP(IDX_1,IDX_2) = 1/MIN_OVERLAP_PIXELS * zscore(bRescaled(fP1:fP1+MIN_OVERLAP_PIXELS-1),1)*zscore(bar2(fP2:fP2+MIN_OVERLAP_PIXELS-1),1)';
                %         pcc

                        [C,IA,IB] = intersect(pos1,pos2);
                        pcc2 = 1/length(IA) * zscore(bRescaled(IA),1)*zscore(bar2(IB),1)';
                        PCC_OVERLAP(IDX_1,IDX_2) = pcc2;
                        len1(IDX_1,IDX_2) = length(b);
                        len2(IDX_1,IDX_2) = length(bar2);
                        lenOverlap(IDX_1,IDX_2) = length(IA);
                else
                end
            % 

        end


    end

end

