function [ rezMax ] = on_compare_rescale_theory2(...
    barcodeGen,theorBar, theorBit, isLinearTF, comparisonMethod,...
    stretchFactors,w,numPixelsAroundBestTheoryMask)
    % on_compare_barcode / in this theorBar opening is moved to the outside
    % loop to speed things up possibly
    % Compares experiments to single theory
    %     Args:
    %         sets: settings structure
    %         barcodeGen: barcode structure
    %         theorBar - theory barcode,theorBit - theory bitmask
    % 
    %     Returns:
    %         comparisonStructure: comparison structure

    % rescales barcode instead of theory
    
    % define a function handle for comparison (i.e. MP or MASSPCC or DTW)
    import CBT.Hca.Core.Comparison.comparison_funs;
    [comparisonFun] = comparison_funs(comparisonMethod,isLinearTF,numPixelsAroundBestTheoryMask);
            
    
    
    rezMaxM = cell(1,length(barcodeGen));
    bestBarStretch = zeros(1,length(barcodeGen));
    bestLength = zeros(1,length(barcodeGen));
    
    rezMax = cell(length(barcodeGen),1);
    % rescale factors are for theory
    
%     [theoryRescale, bitRescale] = nmbp_rescale_direct(theorBar, theorBit, theoryStruct{1}.pixelWidth_nm, theoryStruct{1}.meanBpExt_nm, stretchFactors(j), theoryStruct{1}.psfSigmaWidth_nm, newPsf);

    
    tic
    for i=1:length(barcodeGen)
        i
        barTested = barcodeGen{i}.rawBarcode;
        barBitmask = barcodeGen{i}.rawBitmask;  
        
        lenBarTested = length(barTested);

        for j=1:length(stretchFactors)
            barC = interp1(barTested, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));
            barB = barBitmask(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));
            try
                [rezMax{i}{j}.maxcoef,rezMax{i}{j}.pos,rezMax{i}{j}.or,rezMax{i}{j}.secondPos,rezMax{i}{j}.lengthMatch,~] = comparisonFun(barC, theorBar, barB,theorBit,w);
            catch
                rezMax{i}{j}.maxcoef = 0;rezMax{i}{j}.pos=0;rezMax{i}{j}.or=0;rezMax{i}{j}.secondPos=0;rezMax{i}{j}.lengthMatch=0;rezMax{i}{j}.dist=0;
            end            
        end
        
   end
   toc
   % all max coefs
%    maxCoefs = cellfun(@(x) cellfun(@(y) y.maxcoef(1),x),rezMax,'un',false);
%    xcorrMax(j) = rezMax{j}.maxcoef(1);
   
%    for idx=1:length(barcodeGen)
%        
%    end
%    
%        % find which stretching parameter had the best score
%     [value,b] = max(xcorrMax);
% 
%     % select the results for this best stretching parameter and output
%     % them. If there were no values computed for this barcode, we don't
%     % save anything.
%     if ~isnan(value)
%         rezMaxM{idx} = rezMax{b};
%         rezMaxM{idx}.rescaleFactor = stretchFactors(b);
%         bestBarStretch(idx) = stretchFactors(b);
%         bestLength(idx) = round(lenBarTested*stretchFactors(b));
%         rezMaxM{idx}.allMax = cellfun(@(x) x.maxcoef(1),rezMax);
%         rezMaxM{idx}.allPos = cellfun(@(x) x.pos(1),rezMax);
%         rezMaxM{idx}.allOr = cellfun(@(x) x.or(1),rezMax);
% 
% %             rezMaxM{idx}.bestBarStretch = bestBarStretch;
%     end
%                 
                
                
end 
            


% 
%             % xcorrMax stores the  maximum coefficients
%             xcorrMax = zeros(1,length(stretchFactors));
%             % rezMaz stores the results for one barcode
%             rezMax = cell(1,length(stretchFactors));
% 
%             try
%             % barTested barcode to be tested
%             catch
%             barTested = barcodeGen{idx}.barcode;
%             end
% 
%                
%             % length of this barcode
%             lenBarTested = length(barTested);
% 
%             % barBitmask - bitmask of this barcode
%             try
%             catch
%             barBitmask = barcodeGen{idx}.bitmask;
%             end
% 
%                 % run the loop for the stretch factors
%                 for j=1:length(stretchFactors)
                    % here interpolate both barcode and bitmask 
%                     barC = interp1(barTested, linspace(1,lenBarTested,lenBarTested*stretchFactors(j)));
%                     barB = barBitmask(round(linspace(1,lenBarTested,lenBarTested*stretchFactors(j))));
%                  
%                 end

            
%         end

       
%    end


% % 
% %% test pccs
% pcc = @(x,y) zscore(x,1)'*zscore(y,1)/length(x);
% bar1 = barTested; % or stretched if stretch factor is not 1
% bar2 = theorBar;
% bit1 = barBitmask;
% idxpos = rezMaxM{1}.secondPos(1);
% pos = rezMaxM{1}.pos(1);
% 
% if ~rezMaxM{1}.or
%     frag1 = flipud(bar1);
% else
%     frag1 = bar1;
% end
% % 
% % idxpos = idxpos+1
% frag2 = bar2(pos:pos+length(bar1)-1);
% % 
% bit1 = find(bit1,1,'first');
% % 
% idxpos = idxpos+bit1-1;
% frag1subseq = frag1(idxpos:idxpos+w-1);
% frag2subseq = frag2(idxpos:idxpos+w-1);
% % 
% pcc(frag1subseq',frag2subseq')
