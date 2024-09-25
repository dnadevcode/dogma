function [data,prev] = create_graph_coverage_map(G,overlapStruct1,idGraph)

% load('C:\Users\Lenovo\git\outer\bargraphex.mat')
% 
% find(cellfun(@(x) ~isempty(x),Ggraphs))
% figure,plot(Ggraphs{962})
% 
% G = Ggraphs{962};

if nargin < 3
    idGraph  = 1;
end
% 
% figure
% plot(G,'Layout','force','ArrowSize',5,'MarkerSize',1)


[bin,binsize] = conncomp(G,'Type','weak')

[a,b] = sort(binsize,'desc');

% idGraph = 2;
% binsize = 1;
idx = binsize(bin) == a(idGraph);
SG = subgraph(G, idx);
% figure;plot(SG)
%%
% figure
% plot(SG,'Layout','force','ArrowSize',5,'MarkerSize',1)

%%posCell = cell(1,length(binsize(1)));
%%
prev = {}; % todo: check lengths in every alignment and take average in the end..
data = {};
for i=1:size(SG.Edges,1)
    %
    % get source and sink coordinates 
    source = str2num(SG.Edges.EndNodes{i,1}); % is this source or sink?
    sink = str2num(SG.Edges.EndNodes{i,2}); % is this source or sink?

    %  check if source was discovered before
    sourceGroup =[];
    for j=1:length(prev)
        [intPos,pp] = intersect(prev{j},[source]);
%         intPos
        if ~isempty(intPos)
            sourceGroup = j; 
            pos2 = pp;           

        end
    end
    
    % check if sink was discovered before
	sinkGroup =[];
    for j =1:length(prev)
        [intPos,pp] = intersect(prev{j},[sink]);
        if ~isempty(intPos)
            sinkGroup = j; 
            pos1 = pp;           
        end
    end 
    if sinkGroup==sourceGroup
        % consistency check
    else
 %
    if isempty(sourceGroup)&&isempty(sinkGroup)
        data{end+1} = [overlapStruct1(sink,source).pA overlapStruct1(sink,source).pA+overlapStruct1(sink,source).lenB-1 1;...
            overlapStruct1(sink,source).pB overlapStruct1(sink,source).pB+overlapStruct1(sink,source).lenA-1 overlapStruct1(sink,source).or==1];
        prev{end+1} = [source sink];
    else
        if ~isempty(sourceGroup)
            % define operation which is adding a new element
            newData =  [overlapStruct1(sink,source).pA overlapStruct1(sink,source).pA+overlapStruct1(sink,source).lenB-1 1;...
            overlapStruct1(sink,source).pB overlapStruct1(sink,source).pB+overlapStruct1(sink,source).lenA-1 overlapStruct1(sink,source).or==1];
      
            
%             data{sourceGroup}
            if  data{sourceGroup}(pos2,3)~= newData(1,3)
                newData(1,3) = ~newData(1,3);
                newData(2,3) = ~newData(2,3);

                newData(2,1:2) = [newData(1,2)-newData(2,2)+1 newData(1,2)-newData(2,1)-newData(1,1)+1];
             end

            newData(1:2,1:2) = newData(1:2,1:2)-newData(1,1)+data{sourceGroup}(pos2,1);
            if isempty(sinkGroup) % need to align to newData sink, i.e. newData(2,:)
                data{sourceGroup} =   [ data{sourceGroup}; newData(1,:)];
                prev{sourceGroup} = [ prev{sourceGroup} sink];
            end
            
               if ~isempty(sinkGroup) % need to align to newData sink, i.e. newData(2,:)
                    if  data{sinkGroup}(pos1,3)~= newData(2,3) % if need to flip
                        data{sinkGroup}(:,3)= ~data{sinkGroup}(:,3);
                        data{sinkGroup}(2:end,1:2) = [ data{sinkGroup}(1,2)- data{sinkGroup}(2:end,2)+1  data{sinkGroup}(1,2)- data{sinkGroup}(2:end,1)- data{sinkGroup}(1,1)+1];
                    end
                    data{sinkGroup}(:,1:2) =  data{sinkGroup}(:,1:2)-data{sinkGroup}(pos1,1)+newData(2,1);
                    data{sourceGroup} =   [ data{sourceGroup};  data{sinkGroup}];
                    prev{sourceGroup} = [ prev{sourceGroup} prev{sinkGroup}];
                    data{sinkGroup} = [];
                    prev{sinkGroup} = [];
               end     

        end
        if ~isempty(sinkGroup)&&isempty(sourceGroup)
            % define operation which is adding a new element
            newData =  [overlapStruct1(sink,source).pA overlapStruct1(sink,source).pA+overlapStruct1(sink,source).lenB-1 1;...
            overlapStruct1(sink,source).pB overlapStruct1(sink,source).pB+overlapStruct1(sink,source).lenA-1 overlapStruct1(sink,source).or==1];
      
            
%             data{sinkGroup}
            if  data{sinkGroup}(pos1,3)~= newData(2,3)
                newData(1,3) = ~newData(1,3);
                newData(2,3) = ~newData(2,3);

                newData(2,1:2) = [newData(1,2)-newData(2,2)+1 newData(1,2)-newData(2,1)-newData(1,1)+1];
             end

                newData(1:2,1:2) = newData(1:2,1:2)-newData(2,1)+data{sinkGroup}(pos1,1);
                data{sinkGroup} =   [ data{sinkGroup}; newData(1,:)];
                prev{sinkGroup} = [ prev{sinkGroup} source];


                

        end
     
%         if ~isempty(sinkGroup)&&~isempty(sourceGroup)
% %            newData
%         end

    end
    end
        
        
    end




% end

%%
locs = find(cellfun(@(x) ~isempty(x),data));
figure;hold on
for i=1:size(data{locs(end)},1)
    plot([data{locs(end)}(i,1) data{locs(end)}(i,2)],[i i],'red-')
end

% ix1 =2;
% ix2 =3;
% figure
% plot([data{end}(ix1,1) data{end}(ix1,2)],[1 1],'red.-')
% hold on
% plot([data{end}(ix2,1) data{end}(ix2,2)],[2 2],'red.-')
% ylim([0 3])

    

end

