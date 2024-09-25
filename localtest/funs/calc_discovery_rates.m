
function  [posDif] = calc_discovery_rates(outConsensus,theoryStruct,synthStr,posStructE,islandEltsE)
    
    % 1) align the recovered islands to theory
    islandsGen = cell(1,length(outConsensus));


    % now positions for new
    posNew = cell2mat(cellfun(@(x) x.pos,posStructE,'UniformOutput',0)');
    orDet = cell2mat(cellfun(@(x) x.or,posStructE,'UniformOutput',0)');
    lengthMatch = cell2mat(cellfun(@(x) x.lengthMatch,posStructE,'UniformOutput',0)');

    truePos = cell2mat(cellfun(@(x) x.pos,synthStr,'UniformOutput',0)');
    trueOr = cell2mat(cellfun(@(x) x.or(1),synthStr,'UniformOutput',0)'); % old orientations

    newIndex = cellfun(@(x) x.barid,posStructE);

    truePos = truePos(newIndex); % bars reshuffled
    trueOr = trueOr(newIndex); % bars reshuffled


    for i=1:length(outConsensus)
        isnanmat = ~isnan(outConsensus{i}(1:end,:));

        % coverage
        coverage = sum(isnanmat,1);
        consensus = nanmean(outConsensus{i}(1:end,:));

%         fP = arrayfun(@(x) find(~isnan(outConsensus{i}(x,:)),1,'first'),1:size(outConsensus{i},1));
        islandsGen{i}.rawBarcode = consensus;
%         consensus(:,coverage<=minCoverage) = nan;

        islandsGen{i}.rawBitmask = ones(1,length(consensus));
        islandsGen{i}.rawBitmask(isnan(islandsGen{i}.rawBarcode)) = 0;

        reshufOr = trueOr(islandEltsE{i});
        orClust = orDet(islandEltsE{i});

%         sumV = sum(abs(reshufOr-orClust));
%         sumNew = sum(abs(reshufOr+orClust-3));
%         if sumNew<sumV
%         % need to reorientate vs theory, which means reorient to match the
%         % left
%              islandsGen{i}.rawBarcode = fliplr( islandsGen{i}.rawBarcode );
%              islandsGen{i}.rawBitmask = fliplr( islandsGen{i}.rawBitmask );
%         end

    end

% sF = 0.9:0.01:1.1;


%         % f=figure,plot(find(~isnan(nanmean(outConsensus(2:end,:)))),outConsensus(2:end,~isnan(nanmean(outConsensus(2:end,:))))','.'); 


    sF = 1;%0.95:0.01:1.05;
    sets.comparisonMethod = 'mass_pcc'; % masked mass in case we have complicated mask

    [comparisonStruct] = compare_to_theoretical(islandsGen,theoryStruct,sF,sets);

    % new positions, for found clusters, along theoryStruct. Some do belong
    % to any clusters. Consistency check: clusters should not have repetitive elements ( if it
    % happens, one of them should have been removed)
    newPos = cell2mat(cellfun(@(x) x.pos(1),comparisonStruct,'UniformOutput',0)');
%     newOr = cell2mat(cellfun(@(x) x.or(1),comparisonStruct,'UniformOutput',0)');



    % positions for given cluster
    posDif  = cell(1,length(islandEltsE));
    for ix = 1:length(islandEltsE)
%         [c,d] = sort(b(islandEltsE{ix}));
        posClust = posNew(islandEltsE{ix});
        lengthMatchClust = lengthMatch(islandEltsE{ix});

 

%         orOld(islandEltsE{ix})
        % true positions
        reshufTrue = truePos(islandEltsE{ix});

        reshufOr = trueOr(islandEltsE{ix});
        orClust = orDet(islandEltsE{ix});

        sumV = sum(abs(reshufOr-orClust));
        sumNew = sum(abs(reshufOr+orClust-3));
        if sumNew<sumV
            posClust = -(posClust+lengthMatchClust);%
        end
        
               % same orientation as for consensus bar
        posClust = posClust-min(posClust)+1;

        % todo: also check the orientation
        posOnThry = posClust+newPos(ix)-1;
    
        posDif{ix} = reshufTrue-posOnThry;
    end
%     [extractedPos posOnThry]
%     i=1;
% %     truePos = pos(newIndex( islandElts{i} ));
%     posNew = cell2mat(cellfun(@(x) x.pos,posStructE,'UniformOutput',0)');
%     i = 5
%     newPos = posNew(islandEltsE{i})
% 
% 
%     newIndexE = cellfun(@(x) x.barid,posStructE);
% 
% 
% 
% 
% 
%     for i=1:length(comparisonStruct)
%         truePos = pos(newIndex(islandElts{i}));
%     end
% islandElts
% 
% minOverlap = 150; % test for best min-overlap
% scorethresh = 0.5; % instead use Stouffer like for HC
% tic
% import Core.calc_overlap_pcc_sort_m;
% [overlapStruct] = calc_overlap_pcc_sort_m([barcodeGen], sF,minOverlap,1);
% toc
    % 2) check the recovered positions
    
    % 3) output

end
