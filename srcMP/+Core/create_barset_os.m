function [barcodeIslands, barcodeIslandsData, badData, dataInfo, barIslands] = create_barset_os(sortedVals,sortedIds,barStruct,overlapStruct,synthStr,pxDifSetting)

    % creates bar-set
    %   Args:
    %       pscores - order by which to construct
    %       NN - number of non nan scores
    %       barcodeGen - barcodes
    %       overlapStruct - overlap structure containing the relevant
    %       information about the overlaps between barcodes
    %  

    % todo: also keep the order of merging.
    
    % function is called like this:
    % import Core.create_barset_os;
    % [barcodeIslands,barcodeIslandsData, badData,badDataInfo,barIslands] = ...
    %     create_barset_os(sortedVals(1:lastTrueOverlap-1),sortedIds(1:lastTrueOverlap-1), barcodeGen',oS);

    if nargin < 6
        pxDifSetting = 50;
    end
    
    NN = length(sortedVals);

    % we try add merge barcodes based on the sortedVals. This is combination of
    % creating bargraph and also alignments of barcodes at each step, taking
    % into account the re-scaling
    
     % barcodeIslandsData: start stop orientation resizefactor
    %%
    barcodeIslands = {};
    barcodeIslandsData = {};
    badData = cell(1,NN);
    dataInfo = cell(1,NN);
    barIslands = cell(1,NN);


    [curSinkVec,curSourceVec] =  arrayfun(@(x) ind2sub(size(overlapStruct),sortedIds(x)),1:NN);
    

  
    for idxPair = 1:NN
        % pair to add to barcodeIslands
%         [curSink, curSource] = ind2sub(size(overlapStruct),sortedIds(idxPair));
        if curSinkVec(idxPair)~=0 
        barIslands{idxPair}.curPair =  [curSinkVec(idxPair) curSourceVec(idxPair)];
        curSource = curSourceVec(idxPair);
        curSink = curSinkVec(idxPair);
        
%         % If we want to plot the data
%         import Core.plot_match_simple;
%         [f] = plot_match_simple(barStruct, overlapStruct,curSink,curSource);
%         fig1 = plot_best_pos([], synthStr([curSink curSource]), [], [], [],10000,1);%1/bpPx*10^6
%         if curSink == 4 && curSource == 35 % CAN CHECK SPECIFIC
%             idxPair
%         end

        [sourceGroup,sinkGroup,pos2,pos1] = find_members(barcodeIslands, curSource, curSink);
        % which version of source/sink relation we should use
        tfvec = bin2dec([num2str(isempty(sourceGroup)) num2str(isempty(sinkGroup))]); 

        cS = overlapStruct(curSink,curSource); % should not matter if this is mp overlap or pcc overlap
%         cS.or = -cS.or/2+3/2;    % in case of mp, its reported as 1 -1 
        if cS.or == -1
            cS.or = 2;
        end
       
        % consistency checks: 1) orientation, 2) position, 3) resize factor 4) local score
        dataInfo{idxPair}.positionalDifference = nan;
        switch tfvec
            case 3
                % items haven't appeared before
                barcodeIslandsData{end+1} = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
                -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];
                barcodeIslandsData{end}(:,1:2) = barcodeIslandsData{end}(:,1:2)-barcodeIslandsData{end}(1,1);               
%                  oSpair = test_plot_new_edge(barStruct,  barcodeIslandsData{end},curSink,curSource);
                barcodeIslands{end+1} = [curSource curSink];            
            case 1 % only source is new, put to a source group
                 [barcodeIslandsData, barcodeIslands, newData] = update_barcode_members_info_source(cS, barcodeIslands,barcodeIslandsData,sourceGroup,pos2,curSink);
            case 2  % only sink is new, put to sing group
                [barcodeIslandsData, barcodeIslands, newData] =  update_barcode_members_info_sink(cS,barcodeIslands,barcodeIslandsData,sinkGroup,pos1,curSource);
            case 0
                % most difficult case
%                 newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
%                 -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];

                [barcodeIslandsData, barcodeIslands, badData{idxPair},dataInfo{idxPair}] =  update_barcode_members_info_source_sink(cS,barcodeIslands,barcodeIslandsData,sinkGroup,sourceGroup,pos1,pos2,curSource,curSink,pxDifSetting);
%                     plot_temp(barcodeIslandsData{sourceGroup}, barcodeIslands{sourceGroup},barStruct )
%                 if isfield(dataInfo{idxPair},'failedfilt')
%                     % current curSource - curSink orientation/alignment is
%                     % not good, and we discard all the remaining matches
%                     % including sink/source pair
%                     eltsToRemove = (curSinkVec==curSource)|(curSinkVec==curSink)|  (curSourceVec==curSource)| (curSourceVec==curSink);
%                     curSinkVec(eltsToRemove) = [0];
%                     curSourceVec(eltsToRemove) = [0];
% 
%                 end

            otherwise

        end
        %                  oSpair = test_plot_new_edge(barStruct, newData,curSink,curSource);

%         import Plot.plot_best_pos;
% 
%         fig1 = plot_best_pos([], islandsPosStruct(barcodeIslandsData,barcodeIslands), [], [], [],cumsum(repmat(5000,1,length(barcodeIslandsData))),1);%1/bpPx*10^6

        barIslands{idxPair}.barcodeIslandsData = barcodeIslandsData;
        barIslands{idxPair}.barcodeIslands = barcodeIslands;
        barIslands{idxPair}.posStruct = islandsPosStruct(barcodeIslandsData,barcodeIslands);
%         
%         tableSynth = synth_struct_to_table(synthStr);
%         tableSynth([curSink curSource],:)

        else

        end
        % skip this pair since one of the elements was giving a bad match
        % before
    end
%%

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
% out
% 
% outfull{1}.pos


%     locs = find(cellfun(@(x) ~isempty(x),barcodeIslandsData));
%     figure;hold on
%     for i=1:size(barcodeIslandsData{locs(end)},1)
%         plot([barcodeIslandsData{locs(end)}(i,1) barcodeIslandsData{locs(end)}(i,2)],[i i],'red-')
%     end
%     source = barcodeIslands{locs(end)}(1);
%     sink = barcodeIslands{locs(end)}(2);
%     bS =  barcodeGen{source}.rawBarcode;
%     bS(~barcodeGen{source}.rawBitmask) = nan;
%     figure, plot(barcodeIslandsData{locs(end)}(1,1):barcodeIslandsData{locs(end)}(1,2), bS)
%     hold on
%     bE =  imresize(barcodeGen{sink}.rawBarcode, 'Scale' , [1 barcodeIslandsData{locs(end)}(2,4)]) ;
%     plot(barcodeIslandsData{locs(end)}(2,1):barcodeIslandsData{locs(end)}(2,2), bE(1:end-1))
    
    % conversion to plottable function:
%     oSpair = [];
%     oSpair(curSink,curSource).pA = -(barcodeIslandsData{1}(2,1)-1);
%     oSpair(curSink,curSource).pB = -(barcodeIslandsData{1}(1,1)-1);
%     oSpair(curSink,curSource).lenA = barcodeIslandsData{1}(2,2)-barcodeIslandsData{1}(2,1)+1;
%     oSpair(curSink,curSource).lenB = barcodeIslandsData{1}(1,2)-barcodeIslandsData{1}(1,1)+1;
%     oSpair(curSink,curSource).or =  barcodeIslandsData{1}(2,3);
%     oSpair(curSink,curSource).bestBarStretch =  barcodeIslandsData{1}(2,4);
%     oSpair(curSink,curSource).overlaplen = nan;
%     oSpair(curSink,curSource).score = nan;
%     import Core.plot_match_pcc;
%     %     (1). = 
%     [f] = plot_match_pcc(barStruct, oSpair,curSink,curSource,barStruct);



%     sink =locs(end
    % % 
    % % % barcodeGen = barcodeGen1;
    % % G = digraph;
    % % G = addnode(G,length(barcodeGen));
    % % G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
    % % 
    % % nodeOrientations = ones(1,length(barcodeGen));
    % %  
    % % Ggraphs =cell(1,NN);
    % % 
    % % edges =[];
    % % 
    % % for idxPair=1:NN
    % %     % [maxScore,maxLoc] = max(normalizedScore);
    % %     % idxPair = 1
    % %     badBars = curStd>thrStd;
    % %     if badBars(idxPair)~=1
    % %         [curSink,curSource] = ind2sub(size(PCC_OVERLAP),sortedIds(idxPair));
    % %     
    % %     %     import Core.plot_match_pcc;
    % %     %     [f] = plot_match_pcc(barStruct, overlapStruct,curSink,curSource);
    % %     
    % %         % negative dist mean sink starts to the left of source. 
    % %         dist = overlapStruct(curSink,curSource).pB;
    % %         orEdge = -2*overlapStruct(curSink,curSource).or+3;
    % %     
    % %         try
    % %             path1 = shortestpath(G,curSource,curSink);
    % %         catch
    % %             path1 =[];
    % %         end
    % %         try
    % %             path2 = shortestpath(G,curSink,curSource);
    % %         catch
    % %             path2 = [];
    % %         end
    % %     
    % %        pathExists = ~isempty(path1)||~isempty(path2);
    % %         if pathExists
    % %             matchingOrs = (nodeOrientations(curSource)==orEdge)&&(nodeOrientations(curSink)==1)||...
    % %                 (nodeOrientations(curSource)==(-1*orEdge))&&(nodeOrientations(curSink)==-1);
    % %     %             matchingOrs
    % %             % https://www.pnas.org/doi/pdf/10.1073/pnas.0604040103
    % %             % elimination of false edges with inconsistent orientation
    % %         else
    % %             matchingOrs = 1;
    % %             nodeOrientations(curSink) = orEdge;  
    % %     
    % %         end
    % %     
    % %     
    % %         if matchingOrs
    % %             if dist > 0
    % %                 G = addedge(G, curSource, curSink, 1/dist);
    % %             else
    % %                 G = addedge(G, curSource, curSink, -1/dist); % should be easy to show which one is source.
    % %             end
    % %             edges = [edges; curSource curSink];
    % %         end
    % %     
    % %     
    % %             %
    % %         Gtemp = G;
    % %         
    % %         listN = 1:length(barcodeGen);
    % %         endNodes = Gtemp.Edges.EndNodes(:);
    % %         curNodes = unique(sort(cellfun(@(x) str2num(x),endNodes)'))';
    % %         listN(curNodes) = [];
    % %     
    % %         Gtemp = rmnode(Gtemp,listN);
    % %         weak_bins = conncomp(Gtemp,'Type','weak');
    % %         numComp(idxPair) = max(weak_bins);
    % %         numNodes(idxPair) = length(curNodes);
    % %     
    % %     
    % %     % plot(Gtemp,'Layout','force','ArrowSize',5,'MarkerSize',1)
    % %     
    % %     %     nexttile,plot(Gtemp,'Layout','force','ArrowSize',10,'MarkerSize',2)
    % %     %   title(strcat([' numNodes = ' num2str(length(curNodes)) ', numComp = ' num2str(numComp(i))]))
    % %       
    % %       Ggraphs{idxPair} = Gtemp;
    % %     end
    % % end
    % % 
    % % nonempty = find(cellfun(@(x) ~isempty(x),Ggraphs));
    % % figure
    % % plot(Ggraphs{nonempty(end)},'Layout','force','ArrowSize',5,'MarkerSize',1)
    % % 
    % % finalgraph = Ggraphs{nonempty(end)};

end


function [sourceGroup,sinkGroup,pos2,pos1] = find_members(barcodeIslands,source,sink)

    %  check if source is contained in one of the previous islands
    sourceGroup =[];
    pos2 =[];
    for j=1:length(barcodeIslands)
        % quicker: create a binary vector for each barcodeIsland that has 1
        % at a position of barcode that is contained in barcode island
        [intPos,pp] = intersect(barcodeIslands{j},[source]);

        if ~isempty(intPos)
            sourceGroup = j; 
            pos2 = pp;           
        end
    end
    
      % check if sink was discovered before
	sinkGroup =[];
    pos1 =[];
    for j =1:length(barcodeIslands)
        [intPos,pp] = intersect(barcodeIslands{j},[sink]);
        if ~isempty(intPos)
            sinkGroup = j; 
            pos1 = pp;           
        end
    end 

end

function [barcodeIslandsData, barcodeIslands,newData ] =  update_barcode_members_info_source(cS,barcodeIslands,barcodeIslandsData,sourceGroup,pos2,curSink)

 % define operation which is adding a new element
     newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
        -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];
    newData(:,1:2) =newData(:,1:2)-newData(1,1);

    %
    
%     newData =  [cS.pA cS.pA+cS.lenB-1 1 1;...
%     cS.pB cS.pB+cS.lenA-1 cS.or==1 cS.bestBarStretch];

    sF = barcodeIslandsData{sourceGroup}(pos2,4)/newData(1,4); %  here newData = 1 usually
    newLens = newData(:,2)-newData(:,1)+1;

    newData(2,1) = round(newData(2,1)*sF);
    newData(:,4) = newData(:,4)*sF;
    newData(:,2) = round((newLens)*sF+newData(:,1)-1);

    % check#1 eliminate false positives with incosistent
    % orientation. Relevant if more than one element in the group
    % we're merging.

    
    if barcodeIslandsData{sourceGroup}(pos2,3) ~= newData(1,3)
        newData(1,3) = 3-newData(1,3); % inversion of both
        newData(2,3) = 3-newData(2,3); 
        % also inversion of detected positions (ignoring scaling)
        newData(2,1:2) = [newData(1,2)-newData(2,2) newData(1,2)-newData(2,1)-newData(1,1)];
     end

        newData(1:2,1:2) = newData(1:2,1:2)-newData(1,1)+barcodeIslandsData{sourceGroup}(pos2,1);
        
        % todo: also update resize factor for all elements

%     if isempty(sinkGroup) % need to align to newData sink, i.e. newData(2,:)
        barcodeIslandsData{sourceGroup} =   [ barcodeIslandsData{sourceGroup}; newData(2,:)];
        barcodeIslands{sourceGroup} = [ barcodeIslands{sourceGroup} curSink];
%     end

end


function [barcodeIslandsData, barcodeIslands,newData ] =  update_barcode_members_info_sink(cS,barcodeIslands,barcodeIslandsData,sinkGroup,pos1,source)

	% new data, where we're adding a new source, it's aligned to the group
	% based on 'sink'
    newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
        -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];
    newData(:,1:2) = newData(:,1:2)-newData(1,1); % first row always starts with 0

%     newData =  [cS.pA cS.pA+cS.lenB-1 1 1;...
%     cS.pB cS.pB+cS.lenA-1 cS.or==1 cS.bestBarStretch];

    % check#1 eliminate false positives with incosistent
    % orientation. Relevant if more than one element in the group
    % we're merging.

    % 1) make sure that barcodes re-sized to same factor
%      newData(2,4),     barcodeIslandsData{sinkGroup}(pos1,4)
    sF = barcodeIslandsData{sinkGroup}(pos1,4)/newData(2,4); % keep old, by multiplying by sFold/sFnew
    newLens = newData(:,2)-newData(:,1)+1;
    newData(2,1) = round(newData(2,1)*sF);
    newData(:,4) = newData(:,4)*sF;
    newData(:,2) = round(newLens*sF+newData(:,1)-1);
    
    if barcodeIslandsData{sinkGroup}(pos1,3) ~= newData(2,3)
        newData(:,3) = 3-newData(:,3); % inversion of both
%         newData(2,3) = 3-newData(2,3); 
        % also inversion of detected positions (ignoring scaling)
        newData(2,1:2) = [newData(1,2)-newData(2,2) newData(1,2)-newData(2,1)-newData(1,1)];
%         newData(2,1) = newData(2,1)+6;
%         newData(1,1:2) = [newData(2,2)-newData(1,2)+1 newData(2,2)-newData(1,1)-newData(2,1)+1];
    end


        newData(1:2,1:2) = newData(1:2,1:2)-newData(2,1)+barcodeIslandsData{sinkGroup}(pos1,1);

%     if isempty(sinkGroup) % need to align to newData sink, i.e. newData(2,:)
        barcodeIslandsData{sinkGroup} =   [ barcodeIslandsData{sinkGroup}; newData(1,:)];
        barcodeIslands{sinkGroup} = [ barcodeIslands{sinkGroup} source];
%     end

    %% test: need to plot updated sF barcodes:
%     newData(2,2)-newData(2,1)+1
%     newData(1,2)-newData(1,1)+1
%     oSpair = test_plot_new_edge(barStruct, newData,curSink,curSource);

end


function [barcodeIslandsData, barcodeIslands, badData, badDataInfo ] =  update_barcode_members_info_source_sink(cS,barcodeIslands,barcodeIslandsData,sinkGroup,sourceGroup, pos1,pos2, source, sink, pxDifSetting)
    % in this case both source and sink are available, so we need to check
    % for consistency of orientation, position, and resizing.
    badData = [];
    badDataInfo.positionalDifference = nan;
        newData = [-cS.pB+1 -cS.pB+1+cS.lenB-1  1 1;...
        -cS.pA+1 -cS.pA+1+cS.lenA-1 cS.or cS.bestBarStretch];

    newData(:,1:2) =newData(:,1:2)-newData(1,1);


    tempPosData = newData;
    % Position check: find two elementes in sinkGroup and sourceGroup that
    % are the same. If source and sink in the same group, then easy
    
    % include scaling
    sF = barcodeIslandsData{sourceGroup}(pos2,4)/tempPosData(1,4); %  here newData = 1 usually
    newLens = tempPosData(:,2)-tempPosData(:,1)+1;

    tempPosData(2:end,1) = round(tempPosData(2:end,1)*sF);
    tempPosData(:,4) = tempPosData(:,4)*sF;
    tempPosData(:,2) = round((newLens)*sF+tempPosData(:,1)-1);

    
    % first align sink to source
    if barcodeIslandsData{sourceGroup}(pos2,3) ~= tempPosData(1,3)
        tempPosData(:,3) = 3-tempPosData(:,3); % inversion of both
%         tempPosData(2,3) = 3-tempPosData(2,3); 
        % also inversion of detected positions
        tempPosData(2,1:2) = [tempPosData(1,2)-tempPosData(2,2) tempPosData(1,2)-tempPosData(2,1)-tempPosData(1,1)];
     end

    tempPosData(1:2,1:2) = tempPosData(1:2,1:2)-tempPosData(1,1)+barcodeIslandsData{sourceGroup}(pos2,1);
    
%      oSpair = test_plot_new_edge(barStruct, tempPosData,source,sink);

    %% at which place do we check orientation consistency?
    

    
    
    % check if islands have overlapping elements

    if sinkGroup==sourceGroup % if same group, we already aligned 
        % check orientation!
        orfac = (2*tempPosData(2,3)-3)*(2*barcodeIslandsData{sinkGroup}(pos1,3)-3);
        if  orfac~=1
            badData = [source sink];
            badDataInfo.failedfilt = 1; % failedfilt. 1 means orientation
            return;
        end
        
        %% TODO: test
        % find the positional difference
        positionalDifference = tempPosData(2,1) -  barcodeIslandsData{sinkGroup}(pos1,1);
        badDataInfo.positionalDifference = positionalDifference;
        if abs(positionalDifference) > pxDifSetting % potentially also remove these elements from bargrouping islands altogether since they are maching ambiguosly
            badData = [source sink];
            badDataInfo.failedfilt = 2; % failed positional similarity
            % otherwise,potentially adjust
            return;
        end
    else
        % in this case source and sink is not the same. We align the sink
        % group to the source // we have source in source group and sink in
        % sink group. we merge these groups. There are no elements yet for
        % us to check the positional difference...
        tempDataSink = barcodeIslandsData{sinkGroup};
        
        sF = tempPosData(2,4)/tempDataSink(pos1,4); % keep old, by multiplying by sFold/sFnew
        newLens = tempDataSink(:,2)-tempDataSink(:,1)+1;
        tempDataSink(2:end,1) = round(tempDataSink(2:end,1)*sF);
        tempDataSink(:,4) = tempDataSink(:,4)*sF;
        tempDataSink(:,2) = round(newLens*sF+tempDataSink(:,1)-1);

        if tempDataSink(pos1,3) ~= tempPosData(2,3)
            %% TEST: more than 2 elements
            tempDataSink(:,3) = 3-tempDataSink(:,3); % inversion of both
%             tempDataSink(2,3) = 3-tempDataSink(2,3); 
            % also inversion of detected positions (ignoring scaling)
            tempDataSink(2:end,1:2) = [tempDataSink(1,2)-tempDataSink(2:end,2) tempDataSink(1,2)-tempDataSink(2:end,1)-tempDataSink(1,1)];
    %         newData(2,1) = newData(2,1)+6;
    %         newData(1,1:2) = [newData(2,2)-newData(1,2)+1 newData(2,2)-newData(1,1)-newData(2,1)+1];
        end

        %% TODO: bugchecking: fix+-1 position difference
        tempDataSink(:,1:2) = tempDataSink(:,1:2)-tempDataSink(pos1,1)+tempPosData(2,1);
        

    
        % check if there are pairs of identical elements
        [C,ia,ib] = intersect(barcodeIslands{sourceGroup},barcodeIslands{sinkGroup});

        if isempty(C)
            %% TODO: test
            barcodeIslandsData{sourceGroup} =   [ barcodeIslandsData{sourceGroup};  tempDataSink];
            barcodeIslands{sourceGroup} = [ barcodeIslands{sourceGroup} barcodeIslands{sinkGroup}];
            barcodeIslandsData{sinkGroup} = [];
            barcodeIslands{sinkGroup} = [];
        else
            % check orientation + position
            %% TODO: write up the case when C is non-empty, i.e. two islands share two elements.
%             C
                    % false positives/ Orientation
        % sinkGroup = sourceGroup, mean two barcodes already are in the
        % bargroup
        % check for consistency or orientation 
%         orNew = (2*tempPosData(1,3)-3)*(2*newData(2,3)-3); % 1 if same orientation, -1 if different
%         orPrev =(2*barcodeIslandsData{sinkGroup}(pos1,3)-3)*(2*barcodeIslandsData{sourceGroup}(pos2,3)-3);
%         if orNew ~= orPrev
%             badData = [badData; source sink];
%             badDataInfo.failedfilt = 1; % failedfilt. 1 means orientation
%             return;
%         end

        end
%                 barcodeIslandsData{sinkGroup}(:,1:2) =  barcodeIslandsData{sinkGroup}(:,1:2)-barcodeIslandsData{sinkGroup}(pos1,1)+newData(2,1);

    end

%     
%  
%     
%     if sinkGroup~=sourceGroup
%         % otherwise 
%         if  barcodeIslandsData{sinkGroup}(pos1,3)~= newData(2,3) % if need to flip
%             barcodeIslandsData{sinkGroup}(:,3)= ~barcodeIslandsData{sinkGroup}(:,3);
%             barcodeIslandsData{sinkGroup}(2:end,1:2) = [ barcodeIslandsData{sinkGroup}(1,2)- barcodeIslandsData{sinkGroup}(2:end,2)+1  barcodeIslandsData{sinkGroup}(1,2)- barcodeIslandsData{sinkGroup}(2:end,1)- barcodeIslandsData{sinkGroup}(1,1)+1];
%         end
% 
%         barcodeIslandsData{sinkGroup}(:,1:2) =  barcodeIslandsData{sinkGroup}(:,1:2)-barcodeIslandsData{sinkGroup}(pos1,1)+newData(2,1);
%         barcodeIslandsData{sourceGroup} =   [ barcodeIslandsData{sourceGroup};  barcodeIslandsData{sinkGroup}];
%         barcodeIslands{sourceGroup} = [ barcodeIslands{sourceGroup} barcodeIslands{sinkGroup}];
%         barcodeIslandsData{sinkGroup} = [];
%         barcodeIslands{sinkGroup} = [];
%     end
end


function G=create_graph(barcodeGen,PCC_OVERLAP,sortedIds,NN)
    % first create barcode graph from all:
    G = digraph;
    G = addnode(G,length(barcodeGen));
    G.Nodes.Name = strsplit(num2str(1:length(barcodeGen)))';
    nodeOrientations = ones(1,length(barcodeGen));
    for ii=1:NN
        [curSink,curSource] = ind2sub(size(PCC_OVERLAP),sortedIds(ii));
        G = addedge(G, curSource, curSink, 1); 
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

function oSpair = test_plot_new_edge(barStruct, newData,curSink,curSource);
    oSpair = [];
    oSpair(curSink,curSource).pA = -(newData(2,1)-1);
    oSpair(curSink,curSource).pB = -(newData(1,1)-1);
    oSpair(curSink,curSource).lenA = newData(2,2)-newData(2,1)+1;
    oSpair(curSink,curSource).lenB = newData(1,2)-newData(1,1)+1;
    oSpair(curSink,curSource).or =  newData(2,3);
    oSpair(curSink,curSource).orSource =  newData(1,3);

    oSpair(curSink,curSource).bestBarStretch = newData(2,4);
    oSpair(curSink,curSource).bestBarStretchSource = newData(1,4);

    oSpair(curSink,curSource).overlaplen = nan;
    oSpair(curSink,curSource).score = nan;
    import Core.plot_match_pcc;
    %     (1). = 
    [f] = plot_match_pcc(barStruct, oSpair,curSink,curSource,barStruct);
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