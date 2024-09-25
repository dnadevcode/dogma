function [rezMax, allCoefs] = hca_compare_distance(barcodeGen,theoryStruct, sets )
    % hca_compare_distance
    % Compares experiments to theory. Helping function with the parfor
    %     Args:
    %         sets: settings structure
    %         barcodeGen: barcode structure
    %         theoryStruct: theory structure
    % 
    %     Returns:
    %         comparisonStruct: comparison structure
        
    disp('Starting comparing exp to theory...')
    tic
    
    comparisonMethod = sets.comparisonMethod;
    w = sets.w;
    numPixelsAroundBestTheoryMask = 20;% hardcoded
            
    rezMax = cell(1,length(theoryStruct));
    allCoefs = zeros(length(barcodeGen{1}.rescaled),length(barcodeGen),length(theoryStruct));
    % for the rest of full details without rezMax:
    % allOrs = zeros(length(barcodeGen{1}.rescaled),length(barcodeGen),length(theoryStruct));
    % allLen = zeros(length(barcodeGen{1}.rescaled),length(barcodeGen),length(theoryStruct));
    % allSF = zeros(length(barcodeGen{1}.rescaled),length(barcodeGen),length(theoryStruct));
    % allPos = zeros(length(barcodeGen{1}.rescaled),length(barcodeGen),length(theoryStruct));
    % allPos2 = zeros(length(barcodeGen{1}.rescaled),length(barcodeGen),length(theoryStruct));

    import CBT.Hca.Core.Comparison.on_compare_sf;
    tic
%     h = waitbar(0,'Comparing experiment(s) vs theory(ies)...');

    numOutputs = 1;
    clear futures
    futures(length(theoryStruct)) = parallel.FevalFuture;
    for idx = 1:numel(theoryStruct)
        futures(idx) = parfeval(@on_compare_sf, numOutputs,barcodeGen, theoryStruct(idx),comparisonMethod, w,numPixelsAroundBestTheoryMask);
    end

%     updateWaitbarFutures = afterEach(futures,@(~) waitbar(mean({futures.State} == "finished"),h),0);
%     toc
%     afterAll(updateWaitbarFutures,@(~) delete(h),0);
    wait(futures)

% %     tic
    for idx = 1:length(theoryStruct)
%         [~,value] = fetchNext(futures);
        rezMax{idx} = futures(idx).fetchOutputs{1};
        allCoefs(:,:,idx) =futures(idx).fetchOutputs{2};
    end
    toc

    % Computing distances for each theory against all experiments. This
%     % loop can be parallelized (parfor)
% tic
%    parfor barNr = 1:length(theoryStruct)
% %         barNr
%         try
%             [rezMax{barNr}] = ...
%                 on_compare_sf(barcodeGen,theoryStruct(barNr),comparisonMethod,w,numPixelsAroundBestTheoryMask);
%         catch
%             disp(['error with ',num2str(barNr)]);
%         end
% 
% 
%     end
% toc
    timePassed = toc;
    disp(strcat(['Experiments were compared to theory in ' num2str(timePassed) ' seconds']));
end

