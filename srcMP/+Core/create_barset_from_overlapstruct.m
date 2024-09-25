function [barsetGen] = create_barset_from_overlapstruct(sortedVals, sortedIds, oS, pxDifSetting, sfRatioSetting, barStruct)
% barcodeIslands, barcodeIslandsData, badData, dataInfo, barIslands

    % creates bar-set from overlap structure
    %   Args:
    %       pscores - order by which to construct
    %       NN - number of non nan scores
    %       barcodeGen - barcodes
    %       overlapStruct - overlap structure containing the relevant
    %       information about the overlaps between barcodes
    %  

    %   Returns:
    %       barsetGen

    % todo: also keep the order of merging.
    
    % function is called like this:
    % import Core.create_barset_os;
    % [barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = ...
    %     create_barset_os(sortedVals(1:lastTrueOverlap-1),sortedIds(1:lastTrueOverlap-1), barcodeGen',oS);

    %     % functions:
    %     find_members - which groups does spirce and Reference belong to
    %     update_barcode_members_info_Reference_temp
    %     update_barcode_members_info_Query_temp
    %     update_Query_Reference_temp
    %     check_consistency_new_pair - check consistency of new edge


    if nargin >=6 % convert to struct form cell. Only used for visualizing individual alignments
        if iscell(barStruct)
            barStruct = cell2struct([cellfun(@(x) double(x.rawBarcode),barStruct,'un',false);...
            cellfun(@(x) x.rawBitmask,barStruct,'un',false)]',{'rawBarcode','rawBitmask'},2);
        end

    end

    if nargin < 4
        pxDifSetting = 50; % how much away from correct position we're allowed to be
    end

    if nargin < 5
        sfRatioSetting = 0.02; % how far from correct sf-rescaling ratio we're allowed to be
    end
    
    % keep only non-nan values
    sortedIds = sortedIds(~isnan(sortedVals));
    sortedVals = sortedVals(~isnan(sortedVals));

    NN = length(sortedVals); %total number of values

    % we try add merge barcodes based on the sortedVals. This is combination of
    % creating bargraph and also alignments of barcodes at each step, taking
    % into account the re-scaling
    
    % barcodeIslandsData: start stop orientation resizefactor
    %%
    barcodeIslands = {};
    barcodeIslandsData = {};
    barIslands = cell(1,NN);

    barsetGen.barcodeIslandsStruct = cell(1,NN);
    barsetGen.barcodeIslands = {};
    barsetGen.barcodeIslandsData = {};

    % vectors of things to map
    [curReferenceVec, curQueryVec] =  arrayfun(@(x) ind2sub(size(oS),sortedIds(x)),1:NN);

%     bargroupIdx = 1:NN; % keeps track of bargroup IDX for each cluster group, so that it's easier accessible when needed
  
    for idxPair = 1:NN
        % pair to add to barcodeIslands
%         if curReferenceVec(idxPair)~=0 
        barIslands{idxPair}.curPair =  [curReferenceVec(idxPair) curQueryVec(idxPair)];
        curQuery = curQueryVec(idxPair);
        curReference = curReferenceVec(idxPair);
        
        [QueryGroup,ReferenceGroup,pos2,pos1] = find_members(barsetGen.barcodeIslands, curQuery, curReference);
        % which version of Query/Reference relation we should use
        tfvec = bin2dec([num2str(isempty(QueryGroup)) num2str(isempty(ReferenceGroup))]); 

        barsetGen.currentData{idxPair}.tfvec = tfvec;
        barsetGen.currentData{idxPair}.ReferenceQuery = [curReference  curQuery];
        barsetGen.currentData{idxPair}.QueryReferenceGroup = [QueryGroup ReferenceGroup];
        barsetGen.currentData{idxPair}.accepted = 0;
        cS = oS(curReference,curQuery); % should not matter if this is mp overlap or pcc overlap

        if cS.or == -1
            cS.or = 2;
        end
       
        % consistency checks: 1) orientation, 2) position, 3) resize factor 4) local score
%         dataInfo{idxPair}.positionalDifference = nan;
        barsetGen.consistencyCheck{idxPair}  = [];
        consistencyCheck = [];

        
        switch tfvec
            case 3 % Add new pair  [curQuery curReference] as a new barcodeIsland. Query always starts at 0
                barsetGen.barcodeIslandsData{end+1} = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
                -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];
                barsetGen.barcodeIslandsData{end}(:,1:2) = barsetGen.barcodeIslandsData{end}(:,1:2)-barsetGen.barcodeIslandsData{end}(1,1);               
                barsetGen.barcodeIslands{end+1}.members = [curQuery curReference];  
                barsetGen.barcodeIslands{end}.joingroups = []; % table to keep track if which pairs to merge
        
                barsetGen.barcodeIslands{end}.edges = [idxPair curQuery curReference 3]; 
                barsetGen.currentData{idxPair}.accepted = 1;
                barsetGen.barcodeIslands{end}.bannedList = zeros(1,size(oS,2)); % elements that should not be mapped to a group since they have inconsistent edges
%                 oSpair = test_plot_new_edge(barStruct,  barsetGen.barcodeIslandsData{end},curReference,curQuery);

            case 1 % Reference is new, put to a Query group
                if ~barsetGen.barcodeIslands{QueryGroup}.bannedList(curReference)
                   
                memberPos = []; % check if there is already a previous pair in joingroups
                numPotGroups = size(barsetGen.barcodeIslands{QueryGroup}.joingroups,1);

                if numPotGroups > 0
                    [memberRow,memberPos] = ind2sub([numPotGroups 2],find( ismember(barsetGen.barcodeIslands{QueryGroup}.joingroups(:,3:4),curReference)));
                end
                
                if isempty(memberPos)
                     barsetGen.barcodeIslands{QueryGroup}.joingroups = [   barsetGen.barcodeIslands{QueryGroup}.joingroups ; QueryGroup QueryGroup curQuery curReference]; % same group
                else

%                     barsetGen.barcodeIslands{QueryGroup}.joingroups(memberRow,1:2)

%                     barsetGen.barcodeIslands{QueryGroup}
                    memberRow = memberRow(1);% double check
                    memberPos = memberPos(1);
                    ReferencePast =  barsetGen.barcodeIslands{QueryGroup}.joingroups(memberRow,4); % past Reference
                    QueryPast =  barsetGen.barcodeIslands{QueryGroup}.joingroups(memberRow,3); % past Query
                 
                    if memberPos==1 % add new member as Query
                        posReference = find(ismember( barsetGen.barcodeIslands{QueryGroup}.members,ReferencePast)); % pos in Query group
                        [barSetTemp] =  update_barcode_members_info_Reference_temp(oS(ReferencePast,QueryPast), barsetGen.barcodeIslands{QueryGroup}, barsetGen.barcodeIslandsData{QueryGroup},posReference,QueryPast);
                    else
                         posQuery = find(ismember( barsetGen.barcodeIslands{QueryGroup}.members,QueryPast)); % pos in Query group
                        [barSetTemp] =  update_barcode_members_info_Query_temp(oS(ReferencePast,QueryPast), barsetGen.barcodeIslands{QueryGroup}, barsetGen.barcodeIslandsData{QueryGroup},posQuery,ReferencePast);
                    end
                    % now check consistency with adding curQuery 
                    pos2 = find(ismember(barSetTemp.members,curQuery)); % pos in Query group
                    pos1 = find(ismember(barSetTemp.members,curReference)); % pos in Reference grop

                    % could we re-adjust here to correct slightly for
                    % positional difference?
                    [consistencyCheck, tempPosData ] =  check_consistency_new_pair(cS, barSetTemp.data, pos1, pos2);

                    % maybe tiledlayout figure with consistency check?
%                         oSpair = test_plot_new_edge(barStruct, tempPosData,curReference,curQuery);
                    
%                         plot_temp_island(barSetTemp.data,barSetTemp.members,barStruct )
%                         plot_temp_island(barsetGen.barcodeIslandsData{QueryGroup},barsetGen.barcodeIslands{QueryGroup}.members,barStruct )

                    % if passed consistency checks, update
                    if consistencyCheck.failedfilt ==0 && abs(consistencyCheck.positionalDifference) <= pxDifSetting  && abs(1-consistencyCheck.sFratio) < sfRatioSetting
                        barsetGen.barcodeIslandsData{QueryGroup} = [barSetTemp.data];
                        barsetGen.barcodeIslands{QueryGroup}.members =  barSetTemp.members;
                        barsetGen.barcodeIslands{QueryGroup}.edges = [  barsetGen.barcodeIslands{QueryGroup}.edges;idxPair  QueryPast ReferencePast 1;idxPair curQuery curReference 1]; 
                        barsetGen.currentData{idxPair}.accepted = 1; 
                    else
                        % add curReference to not allowed group so that
                        % it's not merged to this group by accident later
                        barsetGen.barcodeIslands{QueryGroup}.bannedList(curReference) = 1;
                    end
                        barsetGen.barcodeIslands{QueryGroup}.joingroups(memberRow,:) = []; % remove from joingroups
               end
        end
            case 2  % New Query
                if ~barsetGen.barcodeIslands{ReferenceGroup}.bannedList(curQuery)
                  
                memberPos = []; % check if there is already a previous pair in joingroups
                numPotGroups = size(barsetGen.barcodeIslands{ReferenceGroup}.joingroups,1);

                if numPotGroups > 0
                    [memberRow,memberPos] = ind2sub([numPotGroups 2],find( ismember(barsetGen.barcodeIslands{ReferenceGroup}.joingroups(:,3:4),curQuery)));
                end
                            
                if isempty(memberPos) % add new joingroups element
                        barsetGen.barcodeIslands{ReferenceGroup}.joingroups = [barsetGen.barcodeIslands{ReferenceGroup}.joingroups; ReferenceGroup ReferenceGroup curQuery curReference]; % same group
                else
                    memberRow = memberRow(1);
                    memberPos = memberPos(1);
                    % Add previoulsy found edge
                    QueryPast =  barsetGen.barcodeIslands{ReferenceGroup}.joingroups(memberRow,3); % past Query
                    ReferencePast =  barsetGen.barcodeIslands{ReferenceGroup}.joingroups(memberRow,4); % past Reference
                 
                    if memberPos == 1 % Query was new in previous edge
                        posReference = find(ismember( barsetGen.barcodeIslands{ReferenceGroup}.members,ReferencePast)); % pos in Query group
                        [barSetTemp] =  update_barcode_members_info_Reference_temp(oS(ReferencePast,QueryPast), barsetGen.barcodeIslands{ReferenceGroup}, barsetGen.barcodeIslandsData{ReferenceGroup},posReference, QueryPast);
                    else % Reference was new in previous edge
                         posQuery = find(ismember( barsetGen.barcodeIslands{ReferenceGroup}.members,QueryPast)); % pos in Query group
                        [barSetTemp] =  update_barcode_members_info_Query_temp(oS(ReferencePast,QueryPast), barsetGen.barcodeIslands{ReferenceGroup}, barsetGen.barcodeIslandsData{ReferenceGroup},posQuery, ReferencePast);
                    end

                    % now check consistency with adding curQuery 
                    pos2 = find(ismember(barSetTemp.members,curQuery)); % pos in Query group
                    pos1 = find(ismember(barSetTemp.members,curReference)); % pos in Reference grop
                    [consistencyCheck, tempPosData ] =  check_consistency_new_pair(cS, barSetTemp.data, pos1, pos2);

%                     oSpair = test_plot_new_edge(barStruct, tempPosData,curReference,curQuery);
%                     plot_temp_island(barSetTemp.data,barSetTemp.members,barStruct )

                    % if passed consistency checks, update
                    if consistencyCheck.failedfilt ==0 && abs(consistencyCheck.positionalDifference) <= pxDifSetting && abs(1-consistencyCheck.sFratio) < sfRatioSetting
                        barsetGen.barcodeIslandsData{ReferenceGroup} = [barSetTemp.data];
                        barsetGen.barcodeIslands{ReferenceGroup}.members =  barSetTemp.members;
                        barsetGen.barcodeIslands{ReferenceGroup}.edges = [  barsetGen.barcodeIslands{ReferenceGroup}.edges;idxPair QueryPast ReferencePast 2 ;idxPair curQuery curReference 2]; 
                        barsetGen.currentData{idxPair}.accepted = 1;
                    else
                        barsetGen.barcodeIslands{ReferenceGroup}.bannedList(curQuery) = 1;
                    end

                    barsetGen.barcodeIslands{ReferenceGroup}.joingroups(memberRow,:) = []; % remove from joingroups
                end
                end
            case 0 % merge two groups
%                 if QueryGroup==3 && ReferenceGroup ==11
%                     barsetGen.barcodeIslands{11}.joingroups
%                     barsetGen.barcodeIslands{3}.joingroups
% 
%                 end
                if ReferenceGroup == QueryGroup % mean that this pair is already included and consistency checked but we can still run a new consistency check
                        [consistencyCheck, tempPosData ] =  check_consistency_new_pair(cS, barsetGen.barcodeIslandsData{QueryGroup}, pos1, pos2);

                        if consistencyCheck.failedfilt == 0 && abs(consistencyCheck.positionalDifference) <= pxDifSetting  && abs(1-consistencyCheck.sFratio) < sfRatioSetting
                            % add edge to graph
%                             if idxPair == 1016
%                                 idxPair
%                             end
                            barsetGen.barcodeIslands{QueryGroup}.edges = [  barsetGen.barcodeIslands{QueryGroup}.edges;idxPair  curQuery curReference 0]; 
                            barsetGen.currentData{idxPair}.accepted  = 1;
                        else
%                             consistencyCheck % might do some tuning HERE
%                             to improve, since this shows that the
%                             barcode islands needs adjustment
                        % maybe remove this pair
                        end
                else
                % 2) Check if there are sufficient elements to merge in
                % this Reference/Query group

                % First check if there are some elements that should not be
                % merged
                bannedMem1 = sum(barsetGen.barcodeIslands{ReferenceGroup}.bannedList(barsetGen.barcodeIslands{QueryGroup}.members));
                bannedMem2 = sum(barsetGen.barcodeIslands{QueryGroup}.bannedList(barsetGen.barcodeIslands{ReferenceGroup}.members));
                if bannedMem1==0 && bannedMem2 ==0

                memberPosReference = [];
                if size(barsetGen.barcodeIslands{ReferenceGroup}.joingroups,1) > 0 % use bargroupIdx so it's easier to rename things
                    memberPosReference = find( ismember(barsetGen.barcodeIslands{ReferenceGroup}.joingroups(:,2),QueryGroup));
                end
                
                if ~isempty(memberPosReference)% check if there are two unique pairs
                    elts = barsetGen.barcodeIslands{ReferenceGroup}.joingroups(memberPosReference,3:4);
                    if size(elts,1)>= 1

                        C = arrayfun(@(x) isempty(intersect(elts(x,:),[curQuery curReference])),1:size(elts,1));
                        firstUniqueRow = find(C==1);
                        if isempty(firstUniqueRow)
                            memberPosReference = [];
                        else
                            memberPosReference = memberPosReference(firstUniqueRow);
                        end
                    end
                end
                if isempty(memberPosReference)
                     barsetGen.barcodeIslands{ReferenceGroup}.joingroups = [ barsetGen.barcodeIslands{ReferenceGroup}.joingroups ; ReferenceGroup QueryGroup curQuery curReference];
                     barsetGen.barcodeIslands{QueryGroup}.joingroups = [ barsetGen.barcodeIslands{QueryGroup}.joingroups ; QueryGroup ReferenceGroup curQuery curReference]; % might already be an element
                else
                    memberPosReference = memberPosReference(1); % there can be multiple possible merge-edges coming from previous merge. Possibly change that the merge appears at that step already
                    memberPosQuery = find( ismember(barsetGen.barcodeIslands{QueryGroup}.joingroups(:,2),ReferenceGroup));

                    % We update Query Reference with Query-Reference pair
                    QueryReferencePair = barsetGen.barcodeIslands{ReferenceGroup}.joingroups(memberPosReference,3:4);

                    % maybe nicer way ? check which one is Query and which
                    % is Reference
                    pos2 = find(ismember(barsetGen.barcodeIslands{QueryGroup}.members,QueryReferencePair(1))); % pos in Query group
                    pos1 = find(ismember(barsetGen.barcodeIslands{ReferenceGroup}.members,QueryReferencePair(2))); % pos in Reference grop
                    if isempty(pos1)|| isempty(pos2) % change Query
                        pos1 = find(ismember(barsetGen.barcodeIslands{QueryGroup}.members,QueryReferencePair(2))); % pos in Query group
                        pos2 = find(ismember(barsetGen.barcodeIslands{ReferenceGroup}.members,QueryReferencePair(1))); % pos in Reference grop
                        QueryGroupPrev = ReferenceGroup;
                        ReferenceGroupPrev = QueryGroup;
                    else
                        QueryGroupPrev = QueryGroup;
                        ReferenceGroupPrev = ReferenceGroup;
                    end


                    barcodeIslandsReferenceData = barsetGen.barcodeIslandsData{ReferenceGroupPrev};
                    barcodeIslandsQueryData = barsetGen.barcodeIslandsData{QueryGroupPrev};
               
                    % trial merge
%                     idxPair
%                     if idxPair == 137
%                      idxPair
%                     end % this is a bit unclear for the user how we udpdate things.
%                     todo: simplify
                    [barSetTemp,tempDataReference] =  update_Query_Reference_temp(oS(QueryReferencePair(2),QueryReferencePair(1)),barsetGen.barcodeIslands{QueryGroupPrev}.members,barsetGen.barcodeIslands{ReferenceGroupPrev}.members, barcodeIslandsQueryData,barcodeIslandsReferenceData, pos1, pos2);

                    % now check consistency with new pair
                    pos2 = find(ismember(barSetTemp.members,curQuery)); % pos in Query group
                    pos1 = find(ismember(barSetTemp.members,curReference)); % pos in Reference grop
                    
                    % trial merge curQuery curReference
                    [consistencyCheck,tempPosData ] =  check_consistency_new_pair(cS, barSetTemp.data, pos1, pos2);

% 
%                        plot_temp_island(barcodeIslandsQueryData,barsetGen.barcodeIslands{QueryGroupPrev}.members,barStruct )
%                         plot_temp_island(barcodeIslandsReferenceData,barsetGen.barcodeIslands{ReferenceGroupPrev}.members,barStruct )
%                         oSpair = test_plot_new_edge(barStruct, tempPosData,curReference,curQuery);

%                         plot_temp_island(barSetTemp.data,barSetTemp.members,barStruct )

                    if consistencyCheck.failedfilt ==0 && abs(consistencyCheck.positionalDifference) <= pxDifSetting  && abs(1-consistencyCheck.sFratio) < sfRatioSetting
                    % update
                        barsetGen.barcodeIslandsData{QueryGroup} =   [barSetTemp.data];
                        barsetGen.barcodeIslands{QueryGroup}.members =  barSetTemp.members;
                        barsetGen.barcodeIslandsData{ReferenceGroup} = [];
                        % TODO: also update other joingroups in case
                        % there were extra ones in the ReferenceGroup
                        barsetGen.barcodeIslands{QueryGroup}.joingroups(memberPosQuery,:) = [];
                        barsetGen.barcodeIslands{ReferenceGroup}.joingroups(memberPosReference,:) = [];
                        % update other joingroups
                        barsetGen.barcodeIslands{ReferenceGroup}.joingroups(:,1) = QueryGroup; %update Query group
                        barsetGen.barcodeIslands{QueryGroup}.joingroups = [ barsetGen.barcodeIslands{QueryGroup}.joingroups;barsetGen.barcodeIslands{ReferenceGroup}.joingroups ];

                        % rename. This can directly give some mergable
                        % edges
                        for jj=1:size(barsetGen.barcodeIslands{ReferenceGroup}.joingroups,1)
                            barsetGen.barcodeIslands{barsetGen.barcodeIslands{ReferenceGroup}.joingroups(jj,2)}.joingroups(...
                                barsetGen.barcodeIslands{barsetGen.barcodeIslands{ReferenceGroup}.joingroups(jj,2)}.joingroups(:,2)==ReferenceGroup,2)=QueryGroup;
                        end

                        barsetGen.barcodeIslands{QueryGroup}.edges = [  barsetGen.barcodeIslands{QueryGroup}.edges; barsetGen.barcodeIslands{ReferenceGroup}.edges ;idxPair QueryReferencePair 0;idxPair curQuery curReference 0]; 
    
                        barsetGen.barcodeIslands{ReferenceGroup}.members = [];
                        barsetGen.barcodeIslands{ReferenceGroup}.edges = [];
    %                     barsetGen.barcodeIslands{ReferenceGroup}.joingroups = [];
                        barsetGen.currentData{idxPair}.accepted = 1;
                    else % >>
                        barsetGen.barcodeIslands{QueryGroup}.joingroups(memberPosQuery,:) = [];
                        barsetGen.barcodeIslands{ReferenceGroup}.joingroups(memberPosReference,:) = [];
                        barsetGen.barcodeIslands{ReferenceGroup}.bannedList(barsetGen.barcodeIslands{QueryGroup}.members) = 1;
                        barsetGen.barcodeIslands{QueryGroup}.bannedList(barsetGen.barcodeIslands{ReferenceGroup}.members) = 1;
                    end
                end
                else
                    barsetGen.barcodeIslands{ReferenceGroup}.bannedList(barsetGen.barcodeIslands{QueryGroup}.members) = 1;
                    barsetGen.barcodeIslands{QueryGroup}.bannedList(barsetGen.barcodeIslands{ReferenceGroup}.members) = 1;
                end


                end

            otherwise % no case
                               
        end
       barsetGen.consistencyCheck{idxPair} = consistencyCheck;

        %                  oSpair = test_plot_new_edge(barStruct, newData,curReference,curQuery);

%         import Plot.plot_best_pos;
% 
%         fig1 = plot_best_pos([], islandsPosStruct(barcodeIslandsData,barcodeIslands), [], [], [],cumsum(repmat(5000,1,length(barcodeIslandsData))),1);%1/bpPx*10^6

        barIslands{idxPair}.barcodeIslandsData = barcodeIslandsData;
        barIslands{idxPair}.barcodeIslands = barcodeIslands;
        barIslands{idxPair}.posStruct = islandsPosStruct(barcodeIslandsData,barcodeIslands);  
    end
%%

%         % If we want to plot the data
%         import Core.plot_match_simple;
%         [f] = plot_match_simple(barStruct, overlapStruct,curReference,curQuery);
%         fig1 = plot_best_pos([], synthStr([curReference curQuery]), [], [], [],10000,1);%1/bpPx*10^6
%         if curReference == 4 && curQuery == 35 % CAN CHECK SPECIFIC
%             idxPair
%         end


% locs = find(cellfun(@(x) ~isempty(x),barcodeIslandsData));
% idd =4;
% figure;hold on
% for i=1:size(barcodeIslandsData{idd},1)
%     plot([barcodeIslandsData{idd}(i,1) barcodeIslandsData{idd}(i,2)],[i i],'red-')
% end
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], synthStr(barcodeIslands{idd}), [], [], [],10000,1);%1/bpPx*10^6

testplot = 0;
if testplot
import Plot.plot_best_pos;
fig1 = plot_best_pos([], synthStr(cell2mat(barcodeIslands)), [], [], [],10000,1);%1/bpPx*10^6

posStruct = islandsPosStruct(barcodeIslandsData,barcodeIslands);
fig1 = plot_best_pos([], posStruct, [], [], [],cumsum(repmat(2000,1,length(barcodeIslandsData))),1);%1/bpPx*10^6


import Plot.plot_best_pos;
fig1 = plot_best_pos([], synthStr(barcodeIslands{end}), [], [], [],10000,1);%1/bpPx*10^6

posStruct = islandsPosStruct({barcodeIslandsData{end}},barcodeIslands);
fig1 = plot_best_pos([], posStruct, [], [], [],10000,1);%1/bpPx*10^6

badIdx = find(cellfun(@(x) ~isempty(x),dataInfo));
ix = 12;
import Plot.plot_best_pos;
fig1 = plot_best_pos([], synthStr(cell2mat(barIslands{ix}.barcodeIslands)), [], [], [],10000,1);%1/bpPx*10^6

% posStruct = islandsPosStruct(barcodeIslandsData,barcodeIslands);
fig1 = plot_best_pos([], barIslands{ix}.posStruct, [], [], [],cumsum(repmat(2000,1,length(barcodeIslandsData))),1);%1/bpPx*10^6
 out = test_score(barIslands,synthStr,idd)

%%
%  [out,outfull] = test_score(barIslands,synthStr,1,1);

 badPos = [];
 tableSynth = synth_struct_to_table(synthStr);

for i=1:length(barIslands)
    i
    for j=1:length(barIslands{i}.barcodeIslandsData)
        [out,outfull] = test_score(barIslands,tableSynth,j,i);

%         if abs(out.posDif) > 0 %+
%             badPos = [badPos; i j];
%         end
        if  abs(out.orDif)+abs(out.lenDif) > 0 %abs(out.posDif)+
            badPos = [badPos; i j];
        else
        if ~isempty(outfull)
            posdif = abs(diff(outfull{1}.pos'));
            posdif = posdif(posdif>0);
            posdif = posdif(posdif~=10000);

            if ~isempty(posdif)
                badPos = [badPos; i j];
            end
        end
        end
    end
end

badPos(1,:)

ixb = 1

 [out,outfull] = test_score(barIslands,tableSynth, badPos(ixb,2), badPos(ixb,1));
outfull{1}.pos
outfull{1}.len
outfull{1}.or
% abs(diff(outfull{1}.pos'))
barIslands{badPos(ixb,1)}.barcodeIslands{badPos(ixb,2)}
barIslands{badPos(ixb,1)}.barcodeIslandsData{badPos(ixb,2)}

import Plot.plot_best_pos;
fig1 = plot_best_pos([], synthStr(barIslands{badPos(ixb,1)}.barcodeIslands{badPos(ixb,2)}), [], [], [],10000,1);%1/bpPx*10^6

end
end


function [QueryGroup,ReferenceGroup,pos2,pos1] = find_members(barcodeIslands,Query,Reference)

    %  check if Query is contained in one of the previous islands
    QueryGroup =[];
    pos2 =[];
    for j=1:length(barcodeIslands)
        % quicker: create a binary vector for each barcodeIsland that has 1
        % at a position of barcode that is contained in barcode island
        [intPos,pp] = intersect(barcodeIslands{j}.members,[Query]);

        if ~isempty(intPos)
            QueryGroup = j; 
            pos2 = pp;           
        end
    end
    
      % check if Reference was discovered before
	ReferenceGroup =[];
    pos1 =[];
    for j =1:length(barcodeIslands)
        [intPos,pp] = intersect(barcodeIslands{j}.members,[Reference]);
        if ~isempty(intPos)
            ReferenceGroup = j; 
            pos1 = pp;           
        end
    end 

end

function [barSetTemp,newData] =  update_barcode_members_info_Reference_temp(cS,barcodeIslandsReference, barsetGenReferenceTemp,pos1,curQuery)
    %
    if cS.or == -1
        cS.or = 2;
    end
       
	% new data, same for both Reference & Query. First row, for Query, will be
    % added to barsetGenReferenceTemp
    newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
        -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];
    newData(:,1:2) = newData(:,1:2)-newData(1,1); % first row always starts with 0

    
    sF = barsetGenReferenceTemp(pos1,4)/newData(2,4); % keep old, by multiplying by sFold/sFnew
    newLens = newData(:,2)-newData(:,1)+1;
    newData(2,1) = round(newData(2,1)*sF);
    newData(:,4) = newData(:,4)*sF;
    newData(:,2) = round(newLens*sF+newData(:,1)-1);
    
    if barsetGenReferenceTemp(pos1,3) ~= newData(2,3)
        newData(:,3) = 3-newData(:,3); % inversion of both
        newData(2,1:2) = [newData(1,2)-newData(2,2) newData(1,2)-newData(2,1)-newData(1,1)];
    end

    newData(1:2,1:2) = newData(1:2,1:2)-newData(2,1)+barsetGenReferenceTemp(pos1,1);


    barSetTemp.data =   [ barsetGenReferenceTemp; newData(1,:)];
    barSetTemp.members = [barcodeIslandsReference.members curQuery];

end


function [barSetTemp] =  update_barcode_members_info_Query_temp(cS,barcodeIslandsQuery, barsetGenQueryTemp,pos2,curReference)
        if cS.or == -1
            cS.or = 2;
        end
       
 % define operation which is adding a new element
    newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
                -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];
    newData(:,1:2) =newData(:,1:2)-newData(1,1);

    newLens = newData(:,2)-newData(:,1)+1;

    % deal with scaling
    sF = barsetGenQueryTemp(pos2,4)/newData(1,4); %  here newData = 1 usually
    
    newData(2,1) = round(newData(2,1)*sF); % todo: skip rounding in these, and only round when plotting/comparing data, since this introduces errors
    newData(:,4) = newData(:,4)*sF;
    newData(:,2) = round((newLens)*sF+newData(:,1)-1);

    if barsetGenQueryTemp(pos2,3) ~= newData(1,3)
        newData(1,3) = 3-newData(1,3); % inversion of both
        newData(2,3) = 3-newData(2,3); 
        newData(2,1:2) = [newData(1,2)-newData(2,2) newData(1,2)-newData(2,1)-newData(1,1)]; % also inversion of detected positions (ignoring scaling)
     end

    newData(1:2,1:2) = newData(1:2,1:2)-newData(1,1)+barsetGenQueryTemp(pos2,1);
        
    % todo: also update resize factor for all elements
    barSetTemp.data =   [ barsetGenQueryTemp; newData(2,:)];
    barSetTemp.members = [barcodeIslandsQuery.members curReference];

end


%
function [consistencyCheck,tempPosData ] =  check_consistency_new_pair(cS,barcodeIslandsQueryData, pos1, pos2)
    % in this case both Query and Reference are available, so we need to check
    % for consistency of orientation, position, and resizing.
%     badData = [];
%     badDataInfo.positionalDifference = nan;

    % TODO: possibly improve this by comparing all the possible matching
    % pairs?

        if cS.or == -1
            cS.or = 2;
        end
       
    consistencyCheck.failedfilt = 0;


    newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
    -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];

    newData(:,1:2) = newData(:,1:2)-newData(1,1);


    tempPosData = newData;
    % Position check: find two elementes in ReferenceGroup and QueryGroup that
    % are the same. If Query and Reference in the same group, then easy
    
    % include scaling
    sF = barcodeIslandsQueryData(pos2,4)/tempPosData(1,4); %  here newData = 1 usually
    newLens = tempPosData(:,2)-tempPosData(:,1)+1;

    tempPosData(2:end,1) = round(tempPosData(2:end,1)*sF);
    tempPosData(:,4) = tempPosData(:,4)*sF;
    tempPosData(:,2) = round((newLens)*sF+tempPosData(:,1)-1);

    
    % first align Reference to Query
    if barcodeIslandsQueryData(pos2,3) ~= tempPosData(1,3)
        tempPosData(:,3) = 3-tempPosData(:,3); % inversion of both
%         tempPosData(2,3) = 3-tempPosData(2,3); 
        % also inversion of detected positions
        tempPosData(2,1:2) = [tempPosData(1,2)-tempPosData(2,2) tempPosData(1,2)-tempPosData(2,1)-tempPosData(1,1)];
     end

    tempPosData(1:2,1:2) = tempPosData(1:2,1:2)-tempPosData(1,1)+barcodeIslandsQueryData(pos2,1);
    
%      oSpair = test_plot_new_edge(barStruct, tempPosData,Query,Reference);

    % check orientation!
    orfac = (2*tempPosData(2,3)-3)*(2*barcodeIslandsQueryData(pos1,3)-3);
    if  orfac~=1
%         badData = [Query Reference];
        consistencyCheck.failedfilt = 1; % failedfilt. 1 means orientation
    end
    
    consistencyCheck.sFratio = barcodeIslandsQueryData(pos1,4)/tempPosData(2,4);
    
    %% TODO: test
    % find the positional difference
    positionalDifference = tempPosData(2,1) -  barcodeIslandsQueryData(pos1,1);
    consistencyCheck.positionalDifference = positionalDifference;

end

% function [] = check_consistency_new_pair

function [barSetTemp,tempDataReference ] =  update_Query_Reference_temp(cS,barcodeIslandsQuery,barcodeIslandsReference,barcodeIslandsQueryData,barcodeIslandsReferenceData,pos1,pos2)
    % updates Query & Reference groups temporary
    if cS.or == -1
    cS.or = 2;
    end
    newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
    -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];

    newData(:,1:2) = newData(:,1:2)-newData(1,1);

    tempPosData = newData;

    % Position check: find two elementes in ReferenceGroup and QueryGroup that
    % are the same. If Query and Reference in the same group, then easy
    
    % include scaling
    sF = barcodeIslandsQueryData(pos2,4)/tempPosData(1,4); %  here newData = 1 usually
    newLens = tempPosData(:,2)-tempPosData(:,1)+1; % important, the reference element has start 0 in order to have correct re-scaling

    tempPosData(2:end,1) = round(tempPosData(2:end,1)*sF);
    tempPosData(:,4) = tempPosData(:,4)*sF;
    tempPosData(:,2) = round((newLens)*sF+tempPosData(:,1)-1);

    % first align Reference to Query
    if barcodeIslandsQueryData(pos2,3) ~= tempPosData(1,3)
        tempPosData(:,3) = 3-tempPosData(:,3); % inversion of both
        % also inversion of detected positions
        tempPosData(2,1:2) = [tempPosData(1,2)-tempPosData(2,2) tempPosData(1,2)-tempPosData(2,1)-tempPosData(1,1)];
     end

    tempPosData(1:2,1:2) = tempPosData(1:2,1:2)-tempPosData(1,1)+barcodeIslandsQueryData(pos2,1);
    
%      oSpair = test_plot_new_edge(barStruct, tempPosData,Query,Reference);

    % in this case Query and Reference is not the same. We align the Reference
    % group to the Query // we have Query in Query group and Reference in
    % Reference group. we merge these groups. There are no elements yet for
    % us to check the positional difference...
    tempDataReference = barcodeIslandsReferenceData;
    
    sF = tempPosData(2,4)/tempDataReference(pos1,4); % keep old, by multiplying by sFold/sFnew
    newLens = tempDataReference(:,2)-tempDataReference(:,1)+1; % did we scale tempDataReference(pos1,:) to be centered at 0 for correct re-scaling?
    tempDataReference(2:end,1) = round(tempDataReference(2:end,1)*sF);
    tempDataReference(:,4) = tempDataReference(:,4)*sF;
    tempDataReference(:,2) = round(newLens*sF+tempDataReference(:,1)-1);

    if tempDataReference(pos1,3) ~= tempPosData(2,3)
        tempDataReference(:,3) = 3-tempDataReference(:,3); % inversion of both
        % also inversion of detected positions (ignoring scaling)
        tempDataReference(2:end,1:2) = [tempDataReference(1,2)-tempDataReference(2:end,2) tempDataReference(1,2)-tempDataReference(2:end,1)-tempDataReference(1,1)];
    end

    %% TODO: bugchecking: fix+-1 position difference
    tempDataReference(:,1:2) = tempDataReference(:,1:2)-tempDataReference(pos1,1)+tempPosData(2,1);
    
    barSetTemp.data = [barcodeIslandsQueryData;  tempDataReference];
    barSetTemp.members = [barcodeIslandsQuery barcodeIslandsReference];
end




function G=create_graph(barcodeGen,PCC_OVERLAP,sortedIds,NN)
    % first create barcode graph from all:
    G = digraph;
    G = addnode(G,length(barcodeGen));
    G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
    nodeOrientations = ones(1,length(barcodeGen));
    for ii=1:NN
        [curReference,curQuery] = ind2sub(size(PCC_OVERLAP),sortedIds(ii));
        G = addedge(G, curQuery, curReference, 1); 
    end
    Gtemp = G;
    % %         
    listN = 1:length(barcodeGen);
    endNodes = Gtemp.Edges.EndNodes(:);
    curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';
    listN(curNodes) = [];
    % %     
    Gtemp = rmnode(Gtemp,listN);
    % %     weak_bins = conncomp(Gtemp,'Type','weak');
    % %         numComp(idxPair) = max(weak_bins);
    % %         numNodes(idxPair) = length(curNodes);
    % %     
    % %     
    figure;
    plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1)
% %     

end

function oSpair = test_plot_new_edge(barStruct, newData, curReference, curQuery)
    oSpair = [];
    oSpair(curReference,curQuery).pA = -(newData(2,1)-1);
    oSpair(curReference,curQuery).pB = -(newData(1,1)-1);
    oSpair(curReference,curQuery).lenA = newData(2,2)-newData(2,1)+1;
    oSpair(curReference,curQuery).lenB = newData(1,2)-newData(1,1)+1;
    oSpair(curReference,curQuery).or =  newData(2,3);
    oSpair(curReference,curQuery).orQuery =  newData(1,3);

    oSpair(curReference,curQuery).bestBarStretch = newData(2,4);
    oSpair(curReference,curQuery).bestBarStretchQuery = newData(1,4);

    oSpair(curReference,curQuery).overlaplen = nan;
    oSpair(curReference,curQuery).score = nan;
    import Core.plot_match_pcc;
    %     (1). = 
    [f] = plot_match_pcc(barStruct, oSpair,curReference,curQuery,barStruct);
end


function posStruct = islandsPosStruct(barcodeIslandsData,barcodeIslands)
    posStruct = cell(1,size(cell2mat(barcodeIslandsData'),1));
    t= 1;
    for i=1:length(barcodeIslandsData)
        for j=1:size(barcodeIslandsData{i},1)
            posStruct{t}.pos = barcodeIslandsData{i}(j,1);
            posStruct{t}.lengthMatch = barcodeIslandsData{i}(j,2)-barcodeIslandsData{i}(j,1)+1;
            posStruct{t}.or = barcodeIslandsData{i}(j,3);
            posStruct{t}.maxcoef = nan;

            posStruct{t}.idx = i;
            t = t+1;
        end
    end


    
end


function tableSynth = synth_struct_to_table(synthStr)
    tableSynth = zeros(length(synthStr),4);
    for i=1:length(synthStr)
        tableSynth(i,:) = [synthStr{i}.pos synthStr{i}.pos+synthStr{i}.lengthMatch-1 synthStr{i}.or synthStr{i}.rf];
    end

end


function [out,outfull] = test_score(barIslands,tableSynth,idd,ix)
% can evaluate accuracy of single bargrouping island
% idd = 1;
% compare synth to detected
posDif = zeros(1,size(barIslands{ix}.barcodeIslands{idd},1));
lenDif = zeros(1,size(barIslands{ix}.barcodeIslands{idd},1));
orDif = zeros(1,size(barIslands{ix}.barcodeIslands{idd},1));
strFDif = zeros(1,size(barIslands{ix}.barcodeIslands{idd},1));
outfull = [];
for ii=1:size(barIslands{ix}.barcodeIslands{idd},1)
    newData = barIslands{ix}.barcodeIslandsData{idd};
    bars = barIslands{ix}.barcodeIslands{idd};
    alignTo = bars(ii);
    trueData = tableSynth(bars,:);
    
    sF = trueData(ii,4)/newData(ii,4); % keep old, by multiplying by sFold/sFnew
    newLens = newData(:,2)-newData(:,1)+1;
    newData(2:end,1) = round(newData(2:end,1)*sF);
    newData(:,4) = newData(:,4)*sF;
    newData(:,2) = round(newLens*sF+newData(:,1)-1);
    
    if trueData(ii,3) ~= newData(ii,3)
        newData(:,3) = 3-newData(:,3); % inversion of both
%         newData(2,3) = 3-newData(2,3); 
        % also inversion of detected positions (ignoring scaling)
        newData(2:end,1:2) = [newData(1,2)-newData(2:end,2) newData(1,2)-newData(2:end,1)-newData(1,1)];
%         newData(2,1) = newData(2,1)+6;
%         newData(1,1:2) = [newData(2,2)-newData(1,2)+1 newData(2,2)-newData(1,1)-newData(2,1)+1];
    end


    newData(:,1:2) = newData(:,1:2)-newData(ii,1)+trueData(ii,1);

    outfull{ii}.pos =[newData(:,1) trueData(:,1)];
    outfull{ii}.len =[newData(:,2)-newData(:,1)+1 trueData(:,2)-trueData(:,1)+1];
    outfull{ii}.or =[newData(:,3) trueData(:,3)];
    outfull{ii}.strF =[newData(:,4) newData(:,4)];

    posDif(ii) = sqrt(sum((newData(:,1)-trueData(:,1)).^2));
    lenDif(ii) = sqrt(sum((newData(:,2)-newData(:,1)+trueData(:,1)-trueData(:,2)).^2));
    orDif(ii) = sqrt(sum((newData(:,3)-trueData(:,3)).^2));
    strFDif(ii) = sqrt(sum((newData(:,4)-trueData(:,4)).^2));

end
    out.posDif = posDif;
    out.lenDif = lenDif;
    out.orDif = orDif;
    out.strFDif = strFDif;

end




function plot_temp_island(newData,island,barStruct )

    %%
        f=figure;nexttile([1,2]);hold on


    for k=1:length(island)
        barReference = imresize(barStruct(island(k)).rawBarcode,'Scale' ,[1 newData(k,4)]);
        barReferenceBit = imresize(barStruct(island(k)).rawBitmask,'Scale' ,[1 newData(k,4)]);

        if newData(k,3)==2
            barReference = fliplr(barReference);
            barReferenceBit = fliplr(barReferenceBit);
        end
        plot(newData(k,1):newData(k,1)+length(barReference)-1,zscore(barReference))
%         infMat(k,[newData(k,1):newData(k,1)+length(barReference)-1]+10000) = zscore(barReference);

    end
    legend(arrayfun(@(x) num2str(x),island,'un',false))
    
end

% function merge_data()
%     % make code cleaner?
%     ReferencePast =  barsetGen.barcodeIslands{QueryGroup}.joingroups(a,3);
%     QueryPast =  barsetGen.barcodeIslands{QueryGroup}.joingroups(a,4);
%     posReference = find(ismember( barsetGen.barcodeIslands{QueryGroup}.members,ReferencePast)); % pos in Query group
%     [barSetTemp] =  update_barcode_members_info_Reference_temp(oS(ReferencePast,QueryPast), barsetGen.barcodeIslands{QueryGroup}, barsetGen.barcodeIslandsData{QueryGroup},posReference,QueryPast);
% 
%     % now check consistency with adding curQuery 
%     pos2 = find(ismember(barSetTemp.members,curQuery)); % pos in Query group
%     pos1 = find(ismember(barSetTemp.members,curReference)); % pos in Reference grop
% 
%     [consistencyCheck, tempPosData ] =  check_consistency_new_pair(cS, barSetTemp.data, pos1, pos2);
% 
% %                         oSpair = test_plot_new_edge(barStruct, tempPosData,curReference,curQuery);
% %                         plot_temp_island(barSetTemp.data,barSetTemp.members,barStruct )
% 
% 
%     % if passed consistency checks, update
%     if consistencyCheck.failedfilt ==0 && abs(consistencyCheck.positionalDifference) <= pxDifSetting
%         barsetGen.barcodeIslandsData{QueryGroup} = [barSetTemp];
%         barsetGen.barcodeIslands{QueryGroup}.members =  barSetTemp.members;
%     end
%         barsetGen.barcodeIslands{QueryGroup}.joingroups(a,:) = []; % remove from joingroups
% 
% 
% end