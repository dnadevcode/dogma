function [groups,curElts] = consensus_groups(allscores,barcodeGen,Z)

allscores(allscores==0) = nan;
[sortedVals, sortedIds] = sort(allscores(:),'desc','MissingPlacement','last');
[curSinkVec,curSourceVec] =  arrayfun(@(x) ind2sub(size(allscores),sortedIds(x)),1:numel(sortedVals));
curVec = [curSinkVec;curSourceVec];

% get the first merging element of each group
groups = [num2cell(1:length(barcodeGen)) cell(1,size(Z,1))];
curElts = zeros(size(Z,1),2);
stG = length(barcodeGen);
st = 1;
for i = 1:size(Z,1)
    groups{stG+i} = [groups{Z(i,1)} groups{Z(i,2)}];
    % we find first element in curVec that has two elements from this group
    for j=st:size(curVec,2)
        s1 = sum(ismember(curVec(:,j),groups{Z(i,1)}))==1;
        s2 = sum(ismember(curVec(:,j),groups{Z(i,2)}))==1;
        if s1&&s2
            if ismember(curVec(1,j),groups{Z(i,1)})   
                curElts(i,:) = curVec(:,j);
            else
                curElts(i,:) = curVec([2 1],j);
            end
            st = st+1;
            break;
        end
    end
end


end

