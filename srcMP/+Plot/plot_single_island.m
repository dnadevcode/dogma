function plot_single_island(newData,island,barStruct )

    f = figure;nexttile([1,2]);hold on


    for k=1:length(island)
        barSink = imresize(barStruct(island(k)).rawBarcode,'Scale' ,[1 newData(k,4)]);
        barSinkBit = imresize(barStruct(island(k)).rawBitmask,'Scale' ,[1 newData(k,4)]);

        if newData(k,3)==2
            barSink = fliplr(barSink);
            barSinkBit = fliplr(barSinkBit);
        end
        plot(newData(k,1):newData(k,1)+length(barSink)-1,zscore(barSink))
    end
        legend(arrayfun(@(x) num2str(x),island,'un',false))

%     nexttile([1,2]);hold on 
% % 
% % 
% %     %  length(interp1(barStruct(k).rawBarcode, linspace(1, length(barStruct(k).rawBarcode), round(length(barStruct(k).rawBarcode)*overlapStruct(k,iy).bestBarStretch))))
% %     barSink = imresize(barStruct(k).rawBarcode,'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
% %     barSinkBit = imresize(barStruct(k).rawBitmask,'Scale' ,[1 overlapStruct(k,iy).bestBarStretch]);
% %     barSink(~barSinkBit) = nan;
% %     try 
% %         barSource = imresize(barStruct(iy).rawBarcode,'Scale' ,[1 overlapStruct(k,iy).bestBarStretchSource]);
% %         barSourcebit = imresize(barStruct(iy).rawBitmask,'Scale' ,[1 overlapStruct(k,iy).bestBarStretchSource]);
% %         barSource( ~barSourcebit) = nan;
% %     catch
% %         barSource = barStruct2(iy).rawBarcode;
% %         barSource( ~barStruct2(iy).rawBitmask) = nan;
% %     end
% %     pA = overlapStruct(k,iy).pA ;
% %     pB = overlapStruct(k,iy).pB ;
% %     h = overlapStruct(k,iy).overlaplen;
% %     orr = overlapStruct(k,iy).or ;
% %     bS =  overlapStruct(k,iy).bestBarStretch;
% %        
% % 
% %     f=figure;nexttile([1,2]);hold on
% % 
% % if orr==2
% %     barSink = fliplr(barSink);
% % end
% % % if also source flipped
% % try
% %    if  overlapStruct(k,iy).orSource == 2 
% %        barSource = fliplr(barSource);
% %    end
% % catch
% %     
% % end
% % 
% % length(barSink)
% % length(barSource)
% % plot(-pA+1:-pA+length(barSink),nanzscore(barSink)+5,'black')
% % plot(-pB+1:-pB+length(barSource),nanzscore(barSource)+7,'red')
% % 
% % l1 = -pA+1:-pA+length(barSink); l2 = -pB+1:-pB+length(barSource);
% % [C,cc1,cc2] = intersect(l1,l2);
% % 
% % %     pos1 = -pA+1:-pA+length(bBar);
% % %     pos2 = -pB+1:-pB+length(aBar);
% % %     bar2 = b2;
% % 
% % %     fP1 = find(pos1==1,1,'first');
% % %     fP2 = find(pos2==1,1,'first');
% % plot(l1(cc1), nanzscore(barSink(cc1)),'black')
% % plot(l2(cc2), nanzscore(barSource(cc2)),'red')
% % nonnanpos = ~isnan(barSink(cc1)).*(~isnan(barSource(cc2)));
% % b1= barSink(cc1);b2 = barSource(cc2);
% % pcc = nanzscore(b1(find(nonnanpos)))*nanzscore(b2(find(nonnanpos)))'/sum(nonnanpos)
% % 
% % %     pcc = 1/h * zscore(bar1(fP1:fP1+h-1),1)*zscore(bar2(fP2:fP2+h-1),1)';
% % 
% % % 
% % % % do full overlap
% % % lpA = length(bBar); lpB = length(aBar);
% % % 
% % % st = min(pA,pB); % which barcode is more to the left
% % % stop = min(lpA-pA+1,lpB-pB+1);
% % % 
% % % aFul = aBar(pB-st+1:pB+stop-1);
% % % bFul = bBar(pA-st+1:pA+stop-1);
% % % 
% % % plot(-st+1:stop-1,nanzscore( bBar(pA-st+1:pA+stop-1))-6,'black')
% % % plot(-st+1:stop-1, nanzscore(  aBar(pB-st+1:pB+stop-1),1)-6,'red')
% % 
% % % zscore(bBar(pA-st+1:pA+stop-1),1)*zscore(aBar(pB-st+1:pB+stop-1),1)'/length( bBar(pA-st+1:pA+stop-1))
% % 
% % 
% %     text(max(l1(cc1))+50,0,strcat(['alignment C= ' num2str( overlapStruct(k,iy).score) ', len = ' num2str(overlapStruct(k,iy).overlaplen)]))
% % %     text(stop+50,-6,strcat('full overlap C= ',num2str( overlapStruct(k,iy).fullscore)))
% % 
% % % [rescaleFactor, orSign]
% %     try
% %         legend({strcat(['$$bar_{' num2str(k) '}$$ sF='  num2str(bS) ' or= '  num2str(orr) ]),strcat([num2str(iy) ' , sF = ' num2str(overlapStruct(k,iy).bestBarStretchSource)])},'location','southeastoutside','Interpreter','latex');     
% %     catch
% %         legend({strcat(['$$bar_{' num2str(k) '}$$ sF='  num2str(bS) ' or= '  num2str(orr) ]),num2str(iy)},'location','southeastoutside','Interpreter','latex');
% %     end
end