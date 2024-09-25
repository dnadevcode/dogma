function [  ] = plot_comparison_vs_theory(comparisonStruct,theoryStruct,barcodeGen,selectedIndices,savetxt)
    % input len1,selectedIndices,theoryGen,comparisonStructure, names, titleT, savetxt 
    % output ??
    if nargin < 5
        savetxt = 0;
    end

    if savetxt == 1
        defaultMatFilename ='saveplotshere';
        [~, matDirpath] = uiputfile('*.txt', 'Save chosen settings data as', defaultMatFilename);
    end
            
    maxcoef = cell2mat(cellfun(@(x) x.maxcoef,comparisonStruct,'UniformOutput',0)');
    pos = cell2mat(cellfun(@(x) x.pos,comparisonStruct,'UniformOutput',0)');
    orientation = cell2mat(cellfun(@(x) x.or,comparisonStruct,'UniformOutput',0)');

    % Now we plot the barcodes:
    switched = 0;
    import CBT.Hca.UI.Helper.load_theory_and_stretch_ex;
    for i=1:length(selectedIndices)
        ii = selectedIndices(i);
        niceName = theoryStruct{comparisonStruct{ii}.idx}.name;
        pl = [strfind(niceName,'NC'),strfind(niceName,'NZ'),1];
        niceName = niceName(max(pl):end);
        pl = [strfind(niceName,'|'), strfind(niceName,' '),length(niceName)+1];
        niceName = niceName(1:(min(pl)-1));
        
        % first load theory and stretched experiment
        [ theorBar,theorBit, expBar, expBit] = load_theory_and_stretch_ex(ii, theoryStruct, comparisonStruct,barcodeGen );
    
        % if theory was shorter than experiment, then the theory and
        % experiment were flipped positions, and we have to position theory
        % on the experiment, rather than experiment on theory
        
        if length(theorBar) < length(expBar)
           barTemp = theorBar;
           bitTemp = theorBit;
           theorBar = expBar;
           %theorBit = expBit;
           expBar = barTemp;
           expBit = bitTemp;
           switched = 1;
        end
        
        % if different orientation, flip the experimental barcode
        if orientation(ii,1) == 2
            expBar = fliplr(expBar);
            expBit = fliplr(expBit);
        end

        % these are the fit position of experiment on theory
        fitPositions = pos(ii,1):pos(ii,1)+length(expBar)-1;
            
        % find if there are indices that loop over the end of theorBar
        a = find(fitPositions == length(theorBar)+1);
        if ~isempty(a)
            fitPositions(a:end) = 1:length(fitPositions(a:end));
        else
            % find if there are indices that cross 0
            a =  find(fitPositions == 0);
            fitPositions(1:a) = (length(theorBar)-a+1):length(theorBar);
            a = a+1;
        end
    
        barFit = theorBar(fitPositions);

        % mean and std of theory where bitmask of experiment is nonzero.
        % if there were some zero's in the theory, then this would not be
        % correct?
        m1 = mean(barFit(logical(expBit)));
        s1= std(barFit(logical(expBit)));

        % mean and std of experiment
        m2 = mean(expBar(logical(expBit)));
        s2= std(expBar(logical(expBit)));
        
        % plot
        figure;
        hold on;
        plot(theorBar(1:end))
        if ~isempty(a)
            plot([fitPositions(1:a-1) nan fitPositions(a:end)], [((expBar(1:a-1)-m2)/s2)*s1+m1 nan ((expBar(a:end)-m2)/s2)*s1+m1])
        else
            plot(fitPositions, ((expBar-m2)/s2) *s1+m1)
        end   
        xlim([min(fitPositions) max(fitPositions)])

        xlabel('pixel nr.','Interpreter','latex')
        ylabel('Rescaled to theoretical intesity','Interpreter','latex')
        
        if switched == 0
            legend({strrep(niceName,'_','\_'), strcat(['$\hat C_' num2str(ii) '$=' num2str(maxcoef(ii,1),'%0.2f')]) },'Interpreter','latex');
        else
            legend({strcat(['$\hat C=$' num2str(maxcoef(ii,1),'%0.2f')]),niceName},'Interpreter','latex');
        end
        
        if savetxt == 1
            if ii > length(barcodeGen)
                name = 'consensus';
            else
                name = barcodeGen{ii}.name;
            end
            if switched == 0
                CBT.Hca.Export.export_plots_txt(ii,fitPositions,expBar,theorBar,name,niceName,matDirpath)
            else
                warning('Theory was shorter than experiment, so the barcodes were switched places');
                CBT.Hca.Export.export_plots_txt(ii,fitPositions,theorBar,expBar,name,niceName,matDirpath)
            end
        end
    end
end

