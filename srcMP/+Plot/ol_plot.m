function [f] = ol_plot(aBar,bBar, subBarA, subBarB, aFul, bFul, pos,overlapStruct,k,iy,possDif,f)

  % https://se.mathworks.com/help/matlab/ref/bar.html
    if nargin < 12
        f=figure;
        tiledlayout(1,2,'TileSpacing','compact');
        ax1=nexttile([1,2]);hold on
    end

    h = zeros(1,4);
    hold on
    h(1) =plot(possDif+[pos.A(1):pos.A(2)],(aBar-mean(aBar,'omitnan'))./(std(aBar,1,'omitnan')),'Color',  [0 0.4470 0.7410]);
    h(2) = plot(possDif+[pos.B(1):pos.B(2)],(bBar-mean(bBar,'omitnan'))./(std(bBar,1,'omitnan')),'Color',[0.8500 0.3250 0.0980]) ;  
%     figure
%     barh(5,[pos.B(1) pos.B(2)],'stacked')
    if isfield(overlapStruct,'h') % in case no local, don't plot the box
	    h(3) = rectangle('Position', [possDif+pos.lA(1) 5 pos.lA(2)-pos.lA(1)+1 0.5]  , 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.9290 0.6940 0.1250]);% 'LineWidth', 1);
        hline = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color',[0.9290 0.6940 0.1250]);
%     
% 
%         plot(possDif+[pos.lA(1) pos.lA(1)],[(subBarA(1)-mean(aBar))./std(aBar,1) 5],'LineStyle','--','Color','black')
%         plot(possDif+[pos.lA(2) pos.lA(2)],[(subBarA(end)-mean(aBar))./std(aBar,1) 5],'LineStyle','--','Color','black')
    end

    % leftover/global score plot
    try % in case there is left and right left-over
        allvals = pos.Af(1):pos.Af(end);
        if isfield(pos,'lA')
            indicesToRemove = ismember(allvals, pos.lA(1):pos.lA(end));
            allvals(indicesToRemove) = [];
            diff_v = diff(allvals);

            % Find the index where the gap occurs
            gap_index = find(diff_v > 1, 1);

            if ~isempty(gap_index)            
                % Split the vector into two parts at the gap
                part{1} = allvals(1:gap_index);
                part{2} = allvals(gap_index+1:end);
            else
                part{1}  = allvals;
            end
        else
            part{1}  = allvals;
        end

        for jj=1:length(part)
 
            h(4) = rectangle('Position', [possDif+ part{jj}(1)  -5 part{jj}(end)-part{jj}(1)+1 0.5]  , 'FaceColor',[0.4940 0.1840 0.5560],'EdgeColor',[0.4940 0.1840 0.5560]);% 'LineWidth', 1);
            hline = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color',[0.4940 0.1840 0.5560]);
            
            plot(possDif+[part{jj}(1) part{jj}(1)],[(aFul(1)-mean(aBar))./std(aBar,1) -5],'LineStyle','--','Color','black')
            plot(possDif+[part{jj}(end) part{jj}(end)],[(aFul(end)-mean(aBar))./std(aBar,1) -5],'LineStyle','--','Color','black')
        end
    catch
    end

    if isfield(overlapStruct,'h') % in case no local, don't plot the box
%         hline = line(NaN,NaN,'LineWidth',2,'LineStyle','-','Color',[0.9290 0.6940 0.1250]);    
        plot(possDif+[pos.lA(1) pos.lA(1)],[(subBarA(1)-mean(aBar))./std(aBar,1) 5],'LineStyle','--','Color','black')
        plot(possDif+[pos.lA(2) pos.lA(2)],[(subBarA(end)-mean(aBar))./std(aBar,1) 5],'LineStyle','--','Color','black')
    end


    details = 0; % for final figures, we want to reduce information on th eplot
    if details
        % previous, including all CC values etc.
        if isfield(overlapStruct,'h')
            lgd1 = legend({strcat(['$$B_{' num2str(k) '}$$ sF='  num2str(overlapStruct.bestBarStretch) ' or= '  num2str(overlapStruct.or) ]),['$$B_{' num2str(iy),'}$$'],...
                ['$C^{loc,', num2str(overlapStruct.h) , '}_{', num2str(k),',',num2str(iy) , '}$ = ', num2str(overlapStruct.score)],['$C^{leftover,', num2str(overlapStruct.partialLength) , '}_{', num2str(k),',',num2str(iy) , '}$ = ', num2str(overlapStruct.partialScore) ]} ,'location','eastoutside','Interpreter','latex');
        else
            lgd1 = legend({strcat(['$$B_{' num2str(k) '}$$ sF='  num2str(overlapStruct.bestBarStretch) ' or= '  num2str(overlapStruct.or) ]),['$$B_{' num2str(iy),'}$$'],...
                ['Global C=',num2str(overlapStruct.fullscore),',w=',num2str(overlapStruct.overlaplen)] },'location','southeastoutside','Interpreter','latex');
        end
    

    else
        lgd1 = legend({strcat(['$$b_{' num2str(k) '}$$']),['$$b_{' num2str(iy),'}$$'],...
            'local','leftover' },'location','eastoutside','Interpreter','latex');
    end
    



%     lgd1.Layout.Tile = 'east';

end

