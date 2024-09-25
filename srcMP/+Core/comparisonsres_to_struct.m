function [rezStruct2] = comparisonsres_to_struct(barcodeGen, comS, N, minL)
    % 
    %
    %   Args:
    %       barGen, origPos, origFlip, origStr

    if nargin < 4
        minL = 200;
    end

    for k=1:length(barcodeGen)
        comS{k}.maxcoef =   comS{k}.maxcoef(1);
        comS{k}.pos =   comS{k}.pos(1);
        comS{k}.or =   comS{k}.or(1);
    end

%    N = thrS{idxRun}.length ; 
%
rezStruct2 = []; % this outputs "best" pairwise comparison between two barcodes (can be calculated using MP), when both barcodes are
% length re-scaled to their original lengths on the theory barcode
for k=1:length(barcodeGen)
    for iy=1:length(barcodeGen)
        if k~=iy
            rezTemp = comS;
            % make sure that synStemp{i}.pos-synStemp{j}.pos close to each
            % other
            
            % if barcodes differ by more than N, move the barcode which
            % starts further away to the left side. Only relevant in case
            % we have circular barcodes. Todo: implement linear here
            if rezTemp{k}.pos>rezTemp{iy}.pos
                idx = k;
            else
                idx = iy;
            end

            if abs(rezTemp{k}.pos-rezTemp{iy}.pos)>abs(abs(rezTemp{k}.pos-rezTemp{iy}.pos)-N)
                rezTemp{idx}.pos = rezTemp{idx}.pos-N;
            end

                 
            rF1 = (1/comS{iy}.bestBarStretch); % rescale A to B length
            rF2 = 1;%1/origStr(j);
            
            % 

            % length with proper length re-scaling where B is not
            % re-scaled and A re-scaled to B
            lpA =  round(rezTemp{k}.lengthMatch*rF1);
            lpB =  length(barcodeGen{iy}.rawBarcode);%round(comS{iy}.lengthUnrescaled);
         
            rezStruct2(k,iy).posB = 1;
            rezStruct2(k,iy).posA = (rezTemp{k}.pos-rezTemp{iy}.pos); % where barcode starts should be correct, independent on re-scaling of second barcode. Unless one of the barcodes in reversed


            if comS{iy}.or(1) == 2 % should have j always with or=0 and then switch, or allow it to be flipped?
                pT  = (rezTemp{iy}.pos+rezTemp{iy}.lengthMatch-rezTemp{k}.pos-rezTemp{k}.lengthMatch)+1; % circular?
                rezTemp{k}.or =  -rezTemp{k}.or+3;
            else
                pT = rezStruct2(k,iy).posA+1; % 0 should mean identical position
            end
%                 synthStr2(i,j).posA = round((pT)*origStr(j)); % start position is adjusted by origStr(j)
%             else
%%
                % depends if pos or neg. ???
            if pT > 0
                rezStruct2(k,iy).posA =round(pT*(1/comS{iy}.bestBarStretch)); % start position is adjusted by origStr(j) to get correct scaling
            else
                rezStruct2(k,iy).posA =round(pT*(1/comS{iy}.bestBarStretch)); % correct by 1 if before was from 0 
            end
            
            % consider the jth barcode as always at 1 (posB)

            pA =  rezStruct2(k,iy).posA ;
            pB =  rezStruct2(k,iy).posB ;
            
            st = max(pA,pB)-pB+1;%-( pA-pB+1); % start along un-rescaled reference
            stop = min(pA+lpA-1,pB+lpB-1)-pB+1;
            
      
            re = barcodeGen{k}.rawBarcode;
            reb = barcodeGen{k}.rawBitmask;
            reb = reb(round(linspace(1, length(re), length(re)*rF1*comS{k}.bestBarStretch))); % bitmask

            re = interp1(re, linspace(1, length(re), length(re)*rF1*comS{k}.bestBarStretch));
            rt = barcodeGen{iy}.rawBarcode;
            rt = interp1(rt, linspace(1, length(rt), round(length(rt)*rF2)));
            rtb = barcodeGen{iy}.rawBitmask;
            rtb = rtb(round(linspace(1, length(rt), round(length(rt)*rF2))));




            a = [ re re];
            reb = [reb reb];

            b =[ rt rt];
            if rezTemp{k}.or==2
                a = fliplr(a);
                reb = fliplr(reb);
            end
%             if synthStr{j}.or==1 % should have j always with or=0 and then switch, or allow it to be flipped?
%                 b = fliplr(b);
%             end
            t = 0;
            aFul = a(st-rezStruct2(k,iy).posA+1+t:stop-rezStruct2(k,iy).posA+1+t);
            rebFull = reb(st-rezStruct2(k,iy).posA+1+t:stop-rezStruct2(k,iy).posA+1+t);


            bFul = b(st:stop);
            rtbFull = rtb(st:stop);

            pcc = zscore(aFul(find(rebFull.*rtbFull)),1)*zscore(bFul(find(rebFull.*rtbFull)),1)'/length(bFul(find(rebFull.*rtbFull)));

            
            % position of overlap
            rezStruct2(k,iy).pA = st-rezStruct2(k,iy).posA+1 ;
            rezStruct2(k,iy).pB = st ;

            % todo: add leftover score here
            rezStruct2(k,iy).fullscore = pcc;

            rezStruct2(k,iy).overlaplen = sum(rebFull);
            if rezStruct2(k,iy).overlaplen >= minL
                rezStruct2(k,iy).score =   rezStruct2(k,iy).fullscore ;
            else
               rezStruct2(k,iy).score = nan;
            end
            rezStruct2(k,iy).bestBarStretch = rF1*comS{k}.bestBarStretch;
            rezStruct2(k,iy).or = rezTemp{k}.or;
            
            rezStruct2(k,iy).lenB = length(barcodeGen{iy}.rawBarcode);
            rezStruct2(k,iy).lenA = length(interp1(barcodeGen{k}.rawBarcode, linspace(1, length(barcodeGen{k}.rawBarcode), length(barcodeGen{k}.rawBarcode)*rezStruct2(k,iy).bestBarStretch)));
            rezStruct2(k,iy).partialScore = nan;
            rezStruct2(k,iy).partialLength = nan;

        else
            rezStruct2(k,iy).fullscore =nan;
            rezStruct2(k,iy).score =nan;

        end   
    end
end
    
%     k=1;theoryStruct=[];
%     theoryStruct{k}.rawBarcode = refBarcode;
%     theoryStruct{k}.rawBitmask = [];
%     theoryStruct{k}.meanBpExt_nm = 0.3;
%     theoryStruct{k}.psfSigmaWidth_nm = 300;
%     theoryStruct{k}.length = length(refBarcode);
%     theoryStruct{k}.isLinearTF = 0;
%     theoryStruct{k}.name = 'Synthetic theory';
% 


% 
% 
% %% TEST
% test = 0;
% if test
%     rezTemp = comS{idxRun};
%     synthStrT = comS{idxRun};
% 
%     if rezTemp{k}.pos>rezTemp{iy}.pos
%         idx = k;
%     else
%         idx = iy;
%     end
%             
%     if abs(rezTemp{k}.pos-rezTemp{iy}.pos)>abs(abs(rezTemp{k}.pos-rezTemp{iy}.pos)-N)
%         rezTemp{idx}.pos = rezTemp{idx}.pos-N;
%     end
% 
%     % we'll rescale only k'th barcode   
%     rF1 = (synthStrT{iy}.rf); % rescale A to B length
%     rF2 = 1;%1/origStr(j);
%     
%             % 
% 
%     % length with proper length re-scaling where B is not
%     % re-scaled and A re-scaled to B
%     lpA =  round(rezTemp{k}.lengthMatch*rF1);
%     lpB =  round(synthStrT{iy}.lengthUnrescaled);
%          
%     synthStr2T(k,iy).posB = 1;
%     synthStr2T(k,iy).posA = (rezTemp{k}.pos-rezTemp{iy}.pos); % where barcode starts should be correct, independent on re-scaling of second barcode. Unless one of the barcodes in reversed
% 
% 
%     if synthStrT{iy}.or == 2 % should have j always with or=0 and then switch, or allow it to be flipped?
%         pT  = (rezTemp{iy}.pos+rezTemp{iy}.lengthMatch-rezTemp{k}.pos-rezTemp{k}.lengthMatch)+1; % circular?
%         rezTemp{k}.or =  -rezTemp{k}.or+3;
%     else
%         pT = synthStr2T(k,iy).posA+1; % 0 should mean identical position
%     end
% %                 synthStr2(i,j).posA = round((pT)*origStr(j)); % start position is adjusted by origStr(j)
% %             else
% %%
%         % depends if pos or neg. ???
%     if pT > 0
%         synthStr2T(k,iy).posA = round(pT*synthStrT{iy}.rf); % start position is adjusted by origStr(j) to get correct scaling
%     else
%         synthStr2T(k,iy).posA = round(pT*synthStrT{iy}.rf); % correct by 1 if before was from 0 
%     end
%             
%             % consider the jth barcode as always at 1 (posB)
% 
%             
%             pA =  synthStr2T(k,iy).posA ;
%             pB =  synthStr2T(k,iy).posB ;
%             
%             st = max(pA,pB)-pB+1;%-( pA-pB+1); % start along un-rescaled reference
%             stop = min(pA+lpA-1,pB+lpB-1)-pB+1;
%             
%       
%             re = barcodeGen{k}.rawBarcode;
%             reb = barcodeGen{k}.rawBitmask;
%             reb = reb(round(linspace(1, length(re), round(length(re)*rF1*synthStrT{k}.bestBarStretch)))); % bitmask
% 
%             re = interp1(re, linspace(1, length(re), round(length(re)*rF1*synthStrT{k}.bestBarStretch)));
%             rt = barcodeGen{iy}.rawBarcode;
%             rt = interp1(rt, linspace(1, length(rt), round(length(rt)*rF2)));
%             rtb = barcodeGen{iy}.rawBitmask;
%             rtb = rtb(round(linspace(1, length(rt), round(length(rt)*rF2))));
% 
% 
% 
% 
%             a = [ re re];
%             reb = [reb reb];
% 
%             b =[ rt rt];
%             if rezTemp{k}.or==2
%                 a = fliplr(a);
%                 reb = fliplr(reb);
%             end
% %             if synthStr{j}.or==1 % should have j always with or=0 and then switch, or allow it to be flipped?
% %                 b = fliplr(b);
% %             end
%             t = 0;
%             aFul = a(st-synthStr2T(k,iy).posA+1+t:stop-synthStr2T(k,iy).posA+1+t);
%             rebFull = reb(st-synthStr2T(k,iy).posA+1+t:stop-synthStr2T(k,iy).posA+1+t);
% 
% 
%             bFul = b(st:stop);
%             rtbFull = rtb(st:stop);
% 
%             pcc = zscore(aFul(find(rebFull.*rtbFull)),1)*zscore(bFul(find(rebFull.*rtbFull)),1)'/length(bFul(find(rebFull.*rtbFull)))
% 
%             
%             % position of overlap
%             synthStr2T(k,iy).pA = st-synthStr2T(k,iy).posA+1 ;
%             synthStr2T(k,iy).pB = st ;
% 
% 
%             figure;
%             plot(aFul(logical(rebFull.*rtbFull)));
%             hold on
%             plot(bFul(logical(rebFull.*rtbFull)))
%             %%
%             % 
% %             lpA = length(a); lpB = length(b);
% % 
% %             st = min(synthStr2(k,iy).pA,synthStr2(k,iy).pB); % which barcode is more to the left
% %             stop = min(lpA-pA+1,lpB-pB+1);
% %             
% %             aFul = aBar(pB-st+1:pB+stop-1);
% %             bFul = bBar(pA-st+1:pA+stop-1);
% end

end