function [resultStruct] = plot_best_bar_mp(ax, barcodeGen, consensusStruct, comparisonStruct,theoryStruct, maxcoef,onlydirectcomp,sets )
    % plot_best_bar
    
    % plots best barcode vs theory in case barcode is always larger than
    % theory
    if nargin < 7
        onlydirectcomp = 0;
    end
    
 
        %% OLD:
    
    import CBT.Hca.UI.Helper.create_full_table;


    len1=length(barcodeGen);
    % barcode orientations
    orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');

    % number and value of the best barcode
    [dd,ii] =max(maxcoef(:,1));
    
    params = comparisonStruct{ii};
    w = sets.w;
    
    % load theory file
    fileID = fopen(theoryStruct{comparisonStruct{ii}.idx}.filename,'r');
    formatSpec = '%f';
    theorBar = fscanf(fileID,formatSpec);
    fclose(fileID);
    
    
    %%
% %        %% PLOT EXPERIMENTAL BARCODE VS THEORY BARCODE/ SIMPLE REDO - options - linear/circular
% %     % two ways:
% %     % 1) theory is fixed (from px=1 to px=N). 
% %     figure,plot(theorBar)
% %     
% %     figure,plot(1:length(expBar),expBar)
% %     hold on
% %     plot(params.secondPos(1):params.secondPos(1)+length(theorBar)-1,theorBar)
% %     
% %     figure,plot(1:length(expBar),expBar)
% %     hold on
% %     plot(params.secondPos(1):params.secondPos(1)+length(theorBar)-1,flipud(theorBar))
% %     
% %     
% %     expBar = barcodeGen{ii}.rawBarcode;
% %     pq = length(expBar)+(params.secondPos(1)-params.pos(1)+1);
% %     b1 = expBar(pq:pq+params.lengthMatch-1);
% % 
% %     % pq - position on expBar
% % %     pq = abs(params.secondPos(1)-params.pos(1))+1;
% % %     b1 = expBar(pq:pq+params.lengthMatch-1);
% %     if params.or(1)==1
% %         pd =  params.secondPos(1)-(w-length(expBar));
% %     else
% %         pd = params.secondPos(1);
% %     end
% %     b2 = theorBar(pd:pd+params.lengthMatch-1);
% %     plot(pd:pd+params.lengthMatch-1,(b2-nanmean( barfragq{1}))/nanstd(barfragq{1},1),'green')
% % 
% %     if params.or(1)==1
% %         pccV =  pcc(fliplr(b1),b2);
% %     else
% %         pccV =  pcc(b1,b2);
% %     end
% %         
% %     % read off the position from mp, mpI:
% %     import mp.mp_profile_stomp_nan_dna
% %     [maxcoef,pos,or,idxpos,r, mp,mpI,mpPos,mpPIQ] = mp_profile_stomp_nan_dna(expBar', theorBar, ones(1,length(expBar)), ones(1,length(theorBar)), 100, 2^10,50);
% %     
% %     [a,b] = max(mp)
% %     length(expBar)-w
% %     
% %     v1=expBar(b:b+w-1);
% %     import Comparison.pcc;
% %     pos2 = mpI(b);
% %     v2 = theorBar(pos2+ length(expBar)-w:pos2+ length(expBar)-1);
% %     pcc(fliplr(v1),v2)
% %     
% %      abs(params.secondPos(1)-params.pos(1))+1
% %     
     %% now normal
    
    niceName = theoryStruct{comparisonStruct{ii}.idx}.name;
    pl = [strfind(niceName,'NC') strfind(niceName,'NZ')];
    niceName = niceName(pl:end);
    pl = [strfind(niceName,'|') strfind(niceName,' ')];
    niceName = strrep(niceName(1:(min(pl)-1)),'_','\_');
    if isempty(niceName)
        niceName = strrep(theoryStruct{comparisonStruct{ii}.idx}.name,'_','\_');
    end
        
    
    % theory length
    thrLen = theoryStruct{comparisonStruct{ii}.idx}.length;
    
    % bitmask. In case of linear barcode, would like to modify this
    theorBit = ones(1,thrLen);
    
    % load either theory barcode or the consensus barcode
    try
        expBar = barcodeGen{ii}.rawBarcode;
        expBit = barcodeGen{ii}.rawBitmask;
    catch
        try
            expBar = consensusStruct.rawBarcode;
            expBit = consensusStruct.rawBitmask;  
        catch
            expBar = consensusStruct{ii-length(barcodeGen)}.rawBarcode;
            expBit = consensusStruct{ii-length(barcodeGen)}.rawBitmask;  
        end
        
    end

    expLen = length(expBar);

    % interpolate to the length which gave best CC value
    expBar = interp1(expBar, linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch));
    expBit = expBit(round(linspace(1,expLen,expLen*comparisonStruct{ii}.bestBarStretch)));
    expBar(~expBit)= nan;
    
    
    

    
    % in this case we change method that we use for plotting
    resultStruct.bar1 = expBar; % this also needs to be rescaled in case there was scaling
    resultStruct.bar2 = [theorBar'];
    resultStruct.pass = 1;
    shift = find(expBit,1,'first');
%     comparisonStruct{ii}.secondPos = comparisonStruct{ii}.secondPos - shift;
    bar2 =  [resultStruct.bar2  resultStruct.bar2];

    import CBT.Hca.Core.Comparison.pcc;
% %     for i=5:-1:-5
%         fixx=shift-1;
%         b1 = fliplr(resultStruct.bar1(comparisonStruct{2}.secondPos+fixx:comparisonStruct{2}.secondPos+sets.w-1+fixx));
%         b2 = bar2(comparisonStruct{2}.pos+comparisonStruct{2}.secondPos+fixx-1:comparisonStruct{2}.pos+comparisonStruct{2}.secondPos+fixx+sets.w-1-1);
%         pcc(b1,b2)
%     end

    %
     resultStruct.matchTable = [ comparisonStruct{ii}.pos(1)  comparisonStruct{ii}.pos(1)+length(resultStruct.bar1)-1 1 length(resultStruct.bar1) comparisonStruct{ii}.or(1)+1];
     if resultStruct.matchTable(5) == 2
         resultStruct.matchTable(3:4) = resultStruct.matchTable(4:-1:3);
     end
     
%     resultStruct.matchTable = [comparisonStruct{ii}.secondPos(1)+shift-1 comparisonStruct{ii}.secondPos(1)+comparisonStruct{ii}.lengthMatch-1+shift-1 ...
%         comparisonStruct{ii}.pos(1)+comparisonStruct{ii}.secondPos(1)-1+shift-1 comparisonStruct{ii}.pos(1)+comparisonStruct{ii}.secondPos(1)-1+shift-1+comparisonStruct{ii}.lengthMatch-1 ~comparisonStruct{ii}.or(1)+1];
%     resultStruct.matchTable(3:4) = mod(resultStruct.matchTable(3:4)-1,length(resultStruct.bar2))+1;
%     resultStruct.matchTable(1:2) = mod(resultStruct.matchTable(1:2)-1,length(resultStruct.bar1))+1;
%     if resultStruct.matchTable(5) == 2
%         tmp = resultStruct.matchTable(3:4);
%         resultStruct.matchTable(3) = tmp(2);
%         resultStruct.matchTable(4) = tmp(1);
%     end

    
    if onlydirectcomp
        
        matchTableO = resultStruct.matchTable;

        % more complicated in case of circular
        nansLeft = 0;
        if resultStruct.matchTable(1) <= 0
            nansLeft = -resultStruct.matchTable(1)+1;
            resultStruct.matchTable(1) = [1];
            if resultStruct.matchTable(5) == 1 
                resultStruct.matchTable(3) = resultStruct.matchTable(3)+nansLeft;   
            else
                resultStruct.matchTable(3) = resultStruct.matchTable(3)-nansLeft;   

        %             resultStruct.matchTable(4) = resultStruct.matchTable(3)+nansLeft;   

            end
        end
    

        [tempTable,barfragq,barfragr] = create_full_table( resultStruct.matchTable,[theorBar; theorBar],expBar,1);

        barfragq{1} = [nan(nansLeft,1); barfragq{1}];
        barfragr{1} = [nan(1,nansLeft) barfragr{1}];
        
        plot( matchTableO(1):matchTableO(2), (barfragq{1}-nanmean( barfragq{1}))/nanstd(barfragq{1},1),'color','red')
        hold on 
        plot( matchTableO(1):matchTableO(2), (barfragr{1}-nanmean( barfragr{1}))/nanstd(barfragr{1},1),'black')
           
        bpPerPx = theoryStruct{1}.pixelWidth_nm/params.bestBarStretch;

        labelstr = 'Position (px)';
%                 labelstr = 'Position (kbps)';
%         str2 = 'kbps';
%         ticks = 1:50:2*length(bar2);
%         ticksx = floor(ticks*bpPerPx/1000);
%         ax.XTick = [ticks];
%         ax.XTickLabel = [ticksx];   
        xlabel( ax,labelstr,'FontSize', 10,'Interpreter','latex')
    
       title(strcat(['Experimental barcode vs theoretical ']),'Interpreter','latex');
%     %
        % here we actually took position based on all the scores, but it
        % should already be known from the data?
        
%         tmpPar = params.secondPos(1)-params.pos(1);
%         if tmpPar <= 0
%             tmpPar = tmpPar+length(expBar);
%         end
%         pq = tmpPar+1;
        pq = abs(params.secondPos(1)-params.pos(1))+1;
        b1 = expBar(pq:pq+params.lengthMatch-1);
        if params.or(1)==1
            pd =  params.secondPos(1)-(w-length(expBar));
        else
            pd = params.secondPos(1);
        end
        b2 = theorBar(pd:pd+params.lengthMatch-1);
        plot(pd:pd+params.lengthMatch-1,(b2-nanmean( barfragq{1}))/nanstd(barfragq{1},1),'green')

        if params.or(1)==1
            pccV =  pcc(fliplr(b1),b2);
        else
            pccV =  pcc(b1,b2);
        end
%         isequal(params.maxcoef(1), pccV)
%         
%         % we should also plot best position..
%         bpPerPx = 1;
%         labelstr = 'Position (px)';
%         str2 = 'px';
%         ticks = 1:50/bpPerPx:2*length(theorBar);
%         ticksx = resultStruct.matchTable(3)+floor(ticks*bpPerPx);
%         ax.XTick = [ticks];
%         ax.XTickLabel = [ticksx];   
%         xlabel( ax,labelstr,'FontSize', 10,'Interpreter','latex')
%         title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
% %     %
%    
%         [a,b] = max(comparisonStruct{ii}.dist);
%         pq = b;
%         b1 = resultStruct.bar1(pq:pq+comparisonStruct{ii}.lengthMatch-1);
%         pd = comparisonStruct{ii}.secondPos(1);
%         b2 = resultStruct.bar2(pd:pd+comparisonStruct{ii}.lengthMatch-1);
%         plot(pd:pd+comparisonStruct{ii}.lengthMatch-1,(b2-nanmean( barfragq{1}))/nanstd(barfragq{1},1),'green')
% 
%         if comparisonStruct{ii}.or(1)==1
%             pccV =  pcc(fliplr(b1),b2');
%         else
%             pccV =  pcc(b1,b2');
%         end
    
     
     
     
        import CBT.Hca.Core.Comparison.pcc;
        
        legend({strcat(['$\hat C_{\rm ' num2str(ii) '}=$' num2str(pccV,'%0.4f')]), niceName},'Interpreter','latex','location','southoutside')
% 
% 

    else

        import CBT.Hca.UI.Helper.plot_sv_full;
        plot_sv_full(resultStruct,sets,num2str(ii),niceName,0)
    end
     
    
        
%     
%     % expBar with expanded cushion
%     expBar = [ repmat(nan,1,userDefinedSeqCushion) expBar repmat(nan,1,userDefinedSeqCushion)];
% 
% %     numSt = min(userDefinedSeqCushion,userDefinedSeqCushion-comparisonStruct{ii}.pos(1));
%     % don't allow theory start to loop over
%     theoryStart = comparisonStruct{ii}.pos(1)-userDefinedSeqCushion;
%     theoryEnd = comparisonStruct{ii}.pos(1)-userDefinedSeqCushion+length(expBar)-1;
%     
%     if theoryStart < 1
%         % in circular case, these should be taken from the end of theorBar
%         theorBar = [ repmat(nan,abs(theoryStart)+1,1); theorBar];
%         theoryEnd = theoryEnd + abs(theoryStart)+1;
%         theoryStart = 1;
%         thrLen = thrLen+abs(theoryStart)+1;
%     end
%     
%     
%     if theoryEnd > thrLen % now split this into linear and nonlinear case..
% %         try
%             theorBar = [ theorBar; theorBar(1:theoryEnd-thrLen)];
% %         catch
% %             theorBar = [ theorBar; theorBar;theorBar];
% %         end
% %         theorBar = [ theorBar; theorBar; theorBar(1:theoryEnd-thrLen)];
%     end
%     
% 
% 	barStruct.bar1 = expBar;
%     barStruct.bar2 = theorBar;
% 
%     if comparisonStruct{ii}.or(1) == 1
%         barStruct.matchTable = [1 length(expBar) theoryStart theoryEnd ...
%             comparisonStruct{ii}.or(1)  ];
%     else
%         barStruct.matchTable = [1 length(expBar) theoryEnd theoryStart ...
%            comparisonStruct{ii}.or(1)  ];
%     end
%     
%     import CBT.Hca.UI.Helper.create_full_table;
%     [temp_table,barfragq, barfragr] = create_full_table(barStruct.matchTable, barStruct.bar1,barStruct.bar2,1);
% 
%     barfragq{1}(~isnan(barfragq{1})) = zscore(barfragq{1}(~isnan(barfragq{1})));
%     barfragr{1}(~isnan(barfragr{1})) = zscore(barfragr{1}(~isnan(barfragr{1})));
% 
% %     figure,
%     plot(barfragq{1})
%     hold on
%     plot(barfragr{1})
%     
%   
% 
%     xlabel('Position along the sequence cushion (px)','Interpreter','latex')
%     ylabel('Z-scored','Interpreter','latex')
%     if ii <= len1
%         name = num2str(ii);
%     else
%         name = 'consensus';
%     end
% 
%     title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%     %
%     legend({strcat(['$\hat C_{\rm ' name '}=$' num2str(dd,'%0.2f')]), niceName},'Interpreter','latex')
% 
% 

%     
%     
%     % flip barcode if it gabe best match when flipped                  
%     if orientation(ii,1) == 2
%         expBar = fliplr(expBar);
%         expBit = fliplr(expBit);
%     end
%     % fit positions
%     fitPositions = pos(ii,1):pos(ii,1)+length(expBar)-1;
%     
%     if min(fitPositions) <= 0 || max(fitPositions) > thrLen
%         %<= 0 should not happen
% %         max = 
% %         plot(repmat(length(theorBar),1,length(-3:3)),-3:3, 'black');
% %         hold on
%         % the match at best position loops around circularly. Shift the
%         theorBar = [theorBar; theorBar(1:end-1) ];
%         theorBit = [theorBit theorBit(1:end-1) ];
% %         warning('Experiment vs. theory not printed');
%     end
%     
%     % position on the theory barcode
%     barFit = theorBar(fitPositions);
%     barBit = theorBit(fitPositions);
% 
% 
%     % mean and std of theory where bitmask of experiment is nonzero.
%     % if there were some zero's in the theory, then this would not be
%     % correct?
%     m1 = mean(barFit(logical(expBit)));
%     s1= std(barFit(logical(expBit)));
% 
%     % mean and std of experiment
%     m2 = mean(expBar(logical(expBit)));
%     s2= std(expBar(logical(expBit)));
% 
%     % check if score is the same
%     % depends if we have zscore 1 or zscore 0
% %     zscore(barFit(logical(expBit)),1)'*zscore(expBar(logical(expBit)),1)'/length(expBar(logical(expBit)))
% %     dd
% 
%     plot(fitPositions, ((expBar-m2)/s2) *s1+m1)
%     hold on
%     plot(fitPositions, barFit)
%     xlim([min(fitPositions) max(fitPositions)])
% 
%     xlabel('Position (px)','Interpreter','latex')
%     ylabel('Rescaled to theoretical intesity','Interpreter','latex')
%     if ii <= len1
%         name = num2str(ii);
%     else
%         name = 'consensus';
%     end
% 
%     title(strcat(['Experimental barcode vs theory ']),'Interpreter','latex');
%      %
%     legend({strcat(['$\hat C_{\rm ' name '}=$' num2str(dd,'%0.2f')]), niceName},'Interpreter','latex')


  
    
end

