function [synthStr, synthStr2, theoryStruct] = synth_to_struct(barcodeGen, refBarcode, origPos, origFlip, origStr,N,minL)
    % 
    %
    %   Args:
    %       barGen, origPos, origFlip, origStr

refBarcode2 = [refBarcode refBarcode(1:2000)];

% plot the fragment positions. Requires hca function
synthStr = cell(1,length(barcodeGen));
for k=1:length(barcodeGen)
    synthStr{k}.idx = 1;
    synthStr{k}.pos = origPos(k);
    synthStr{k}.lengthMatch = round(length(barcodeGen{k}.rawBarcode)/origStr(k));
    synthStr{k}.lengthUnrescaled = length(barcodeGen{k}.rawBarcode);
    synthStr{k}.bestBarStretch = 1./origStr(k);
    synthStr{k}.or = double(origFlip(k))+1;
    synthStr{k}.rf = origStr(k);

    % calc coef to theory (at best extension)
    re = barcodeGen{k}.rawBarcode;
    reb = barcodeGen{k}.rawBitmask(round(linspace(1, length(re), round(length(re)*synthStr{k}.bestBarStretch))));

    re = interp1(re, linspace(1, length(re), round(length(re)*synthStr{k}.bestBarStretch)));
    
    if  synthStr{k}.or == 2
        re = fliplr(re);
        reb = fliplr(reb);
    end
%     synthStr{i}.maxcoef = 0.99; % calculate max coef against thry

    refBar = refBarcode2(synthStr{k}.pos:synthStr{k}.pos+synthStr{k}.lengthMatch-1);
    synthStr{k}.maxcoef = zscore(refBar(reb),1)*zscore(re(reb),1)'/sum(reb);

%     figure,plot(refBarcode2(synthStr{i}.pos:synthStr{i}.pos+synthStr{i}.lengthMatch-1))
%     hold on
%     plot(re)
end


% Need to redo

% if circ, check if loop around
% synthStrLin = synthStr;
% % for i=1:length(barcodeGen)
% %     if synthStr{i}.pos+synthStr{i}.lengthMatch-1>TOTAL_RAND_LENGTH
% %         synthStrLin{i}.pos = synthStrLin{i}.pos-TOTAL_RAND_LENGTH;
% %     end
% % end
% % expected overlap positions
%
synthStr2 = []; % this outputs "best" pairwise comparison between two barcodes (can be calculated using MP), when both barcodes are
% length re-scaled to their original lengths on the theory barcode
for k=1:length(barcodeGen)
    for iy=1:length(barcodeGen)
        if k~=iy
            synStemp = synthStr;
            % make sure that synStemp{i}.pos-synStemp{j}.pos close to each
            % other
            
            % if barcodes differ by more than N, move the barcode which
            % starts further away to the left side. Only relevant in case
            % we have circular barcodes. Todo: implement linear here
            if synStemp{k}.pos>synStemp{iy}.pos
                idx = k;
            else
                idx = iy;
            end

            if abs(synStemp{k}.pos-synStemp{iy}.pos)>abs(abs(synStemp{k}.pos-synStemp{iy}.pos)-N)
                synStemp{idx}.pos = synStemp{idx}.pos-N;
            end

                 
            rF1 = (synthStr{iy}.rf); % rescale A to B length
            rF2 = 1;%1/origStr(j);
            
            % 

            % length with proper length re-scaling where B is not
            % re-scaled and A re-scaled to B
            lpA =  round(synStemp{k}.lengthMatch*rF1);
            lpB =  round(synthStr{iy}.lengthUnrescaled);
         
            synthStr2(k,iy).posB = 1;
            synthStr2(k,iy).posA = (synStemp{k}.pos-synStemp{iy}.pos); % where barcode starts should be correct, independent on re-scaling of second barcode. Unless one of the barcodes in reversed


            if synthStr{iy}.or == 2 % should have j always with or=0 and then switch, or allow it to be flipped?
                pT  = (synStemp{iy}.pos+synStemp{iy}.lengthMatch-synStemp{k}.pos-synStemp{k}.lengthMatch)+1; % circular?
                synStemp{k}.or =  -synStemp{k}.or+3;
            else
                pT = synthStr2(k,iy).posA+1; % 0 should mean identical position
            end
%                 synthStr2(i,j).posA = round((pT)*origStr(j)); % start position is adjusted by origStr(j)
%             else
%%
                % depends if pos or neg. ???
            if pT > 0
                synthStr2(k,iy).posA =round(pT*origStr(iy)); % start position is adjusted by origStr(j) to get correct scaling
            else
                synthStr2(k,iy).posA =round(pT*origStr(iy)); % correct by 1 if before was from 0 
            end
            
            % consider the jth barcode as always at 1 (posB)

            pA =  synthStr2(k,iy).posA ;
            pB =  synthStr2(k,iy).posB ;
            
            st = max(pA,pB)-pB+1;%-( pA-pB+1); % start along un-rescaled reference
            stop = min(pA+lpA-1,pB+lpB-1)-pB+1;
            
      
            re = barcodeGen{k}.rawBarcode;
            reb = barcodeGen{k}.rawBitmask;
            reb = reb(round(linspace(1, length(re), length(re)*rF1*synthStr{k}.bestBarStretch))); % bitmask

            re = interp1(re, linspace(1, length(re), length(re)*rF1*synthStr{k}.bestBarStretch));
            rt = barcodeGen{iy}.rawBarcode;
            rt = interp1(rt, linspace(1, length(rt), round(length(rt)*rF2)));
            rtb = barcodeGen{iy}.rawBitmask;
            rtb = rtb(round(linspace(1, length(rt), round(length(rt)*rF2))));




            a = [ re re];
            reb = [reb reb];

            b =[ rt rt];
            if synStemp{k}.or==2
                a = fliplr(a);
                reb = fliplr(reb);
            end
%             if synthStr{j}.or==1 % should have j always with or=0 and then switch, or allow it to be flipped?
%                 b = fliplr(b);
%             end
            t = 0;
            aFul = a(st-synthStr2(k,iy).posA+1+t:stop-synthStr2(k,iy).posA+1+t);
            rebFull = reb(st-synthStr2(k,iy).posA+1+t:stop-synthStr2(k,iy).posA+1+t);


            bFul = b(st:stop);
            rtbFull = rtb(st:stop);

            pcc = zscore(aFul(find(rebFull.*rtbFull)),1)*zscore(bFul(find(rebFull.*rtbFull)),1)'/length(bFul(find(rebFull.*rtbFull)));

            
            % position of overlap
            synthStr2(k,iy).pA = st-synthStr2(k,iy).posA+1 ;
            synthStr2(k,iy).pB = st ;


            synthStr2(k,iy).fullscore = pcc;

            synthStr2(k,iy).overlaplen = sum(rebFull);
            if synthStr2(k,iy).overlaplen >= minL
                synthStr2(k,iy).score =   synthStr2(k,iy).fullscore ;
            else
               synthStr2(k,iy).score = nan;
            end
            synthStr2(k,iy).bestBarStretch = rF1*synthStr{k}.bestBarStretch;
            synthStr2(k,iy).or = synStemp{k}.or;
            
            synthStr2(k,iy).lenB = length(barcodeGen{iy}.rawBarcode);
            synthStr2(k,iy).lenA = length(interp1(barcodeGen{k}.rawBarcode, linspace(1, length(barcodeGen{k}.rawBarcode), length(barcodeGen{k}.rawBarcode)*synthStr2(k,iy).bestBarStretch)));
%             synthStr2(i,j).fullscore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
%             overlapStruct(k,iy).overlaplen = length(aFul);
%             overlapPos = 
            synthStr2(k,iy).partialScore = nan;
            synthStr2(k,iy).partialLength = nan;

        else
            synthStr2(k,iy).fullscore =nan;
            synthStr2(k,iy).score =nan;
            synthStr2(k,iy).partialScore = nan;
            synthStr2(k,iy).partialLength = nan;
            synthStr2(k,iy).lenB = nan;
            synthStr2(k,iy).lenA = nan;            
            synthStr2(k,iy).overlaplen = nan;
%             synthStr2(k,iy).partialLength = nan;

        end   
    end
end
    
    k=1;theoryStruct=[];
    theoryStruct{k}.rawBarcode = refBarcode;
    theoryStruct{k}.rawBitmask = [];
    theoryStruct{k}.meanBpExt_nm = 0.3;
    theoryStruct{k}.psfSigmaWidth_nm = 300;
    theoryStruct{k}.length = length(refBarcode);
    theoryStruct{k}.isLinearTF = 0;
    theoryStruct{k}.name = 'Synthetic theory';





%% TEST
test = 0;
if test
    synStemp = synthStr{idxRun};
    synthStrT = synthStr{idxRun};

    if synStemp{k}.pos>synStemp{iy}.pos
        idx = k;
    else
        idx = iy;
    end
            
    if abs(synStemp{k}.pos-synStemp{iy}.pos)>abs(abs(synStemp{k}.pos-synStemp{iy}.pos)-N)
        synStemp{idx}.pos = synStemp{idx}.pos-N;
    end

    % we'll rescale only k'th barcode   
    rF1 = (synthStrT{iy}.rf); % rescale A to B length
    rF2 = 1;%1/origStr(j);
    
            % 

    % length with proper length re-scaling where B is not
    % re-scaled and A re-scaled to B
    lpA =  round(synStemp{k}.lengthMatch*rF1);
    lpB =  round(synthStrT{iy}.lengthUnrescaled);
         
    synthStr2T(k,iy).posB = 1;
    synthStr2T(k,iy).posA = (synStemp{k}.pos-synStemp{iy}.pos); % where barcode starts should be correct, independent on re-scaling of second barcode. Unless one of the barcodes in reversed


    if synthStrT{iy}.or == 2 % should have j always with or=0 and then switch, or allow it to be flipped?
        pT  = (synStemp{iy}.pos+synStemp{iy}.lengthMatch-synStemp{k}.pos-synStemp{k}.lengthMatch)+1; % circular?
        synStemp{k}.or =  -synStemp{k}.or+3;
    else
        pT = synthStr2T(k,iy).posA+1; % 0 should mean identical position
    end
%                 synthStr2(i,j).posA = round((pT)*origStr(j)); % start position is adjusted by origStr(j)
%             else
%%
        % depends if pos or neg. ???
    if pT > 0
        synthStr2T(k,iy).posA = round(pT*synthStrT{iy}.rf); % start position is adjusted by origStr(j) to get correct scaling
    else
        synthStr2T(k,iy).posA = round(pT*synthStrT{iy}.rf); % correct by 1 if before was from 0 
    end
            
            % consider the jth barcode as always at 1 (posB)

            
            pA =  synthStr2T(k,iy).posA ;
            pB =  synthStr2T(k,iy).posB ;
            
            st = max(pA,pB)-pB+1;%-( pA-pB+1); % start along un-rescaled reference
            stop = min(pA+lpA-1,pB+lpB-1)-pB+1;
            
      
            re = barcodeGen{k}.rawBarcode;
            reb = barcodeGen{k}.rawBitmask;
            reb = reb(round(linspace(1, length(re), round(length(re)*rF1*synthStrT{k}.bestBarStretch)))); % bitmask

            re = interp1(re, linspace(1, length(re), round(length(re)*rF1*synthStrT{k}.bestBarStretch)));
            rt = barcodeGen{iy}.rawBarcode;
            rt = interp1(rt, linspace(1, length(rt), round(length(rt)*rF2)));
            rtb = barcodeGen{iy}.rawBitmask;
            rtb = rtb(round(linspace(1, length(rt), round(length(rt)*rF2))));




            a = [ re re];
            reb = [reb reb];

            b =[ rt rt];
            if synStemp{k}.or==2
                a = fliplr(a);
                reb = fliplr(reb);
            end
%             if synthStr{j}.or==1 % should have j always with or=0 and then switch, or allow it to be flipped?
%                 b = fliplr(b);
%             end
            t = 0;
            aFul = a(st-synthStr2T(k,iy).posA+1+t:stop-synthStr2T(k,iy).posA+1+t);
            rebFull = reb(st-synthStr2T(k,iy).posA+1+t:stop-synthStr2T(k,iy).posA+1+t);


            bFul = b(st:stop);
            rtbFull = rtb(st:stop);

            pcc = zscore(aFul(find(rebFull.*rtbFull)),1)*zscore(bFul(find(rebFull.*rtbFull)),1)'/length(bFul(find(rebFull.*rtbFull)))

            
            % position of overlap
            synthStr2T(k,iy).pA = st-synthStr2T(k,iy).posA+1 ;
            synthStr2T(k,iy).pB = st ;


            figure;
            plot(aFul(logical(rebFull.*rtbFull)));
            hold on
            plot(bFul(logical(rebFull.*rtbFull)))
            %%
            % 
%             lpA = length(a); lpB = length(b);
% 
%             st = min(synthStr2(k,iy).pA,synthStr2(k,iy).pB); % which barcode is more to the left
%             stop = min(lpA-pA+1,lpB-pB+1);
%             
%             aFul = aBar(pB-st+1:pB+stop-1);
%             bFul = bBar(pA-st+1:pA+stop-1);
end

end