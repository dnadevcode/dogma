
function [barcodeGen, synthStr, synthStr2, theoryStruct, refBarcode, origPos, origFlip,origStr] = ...
    gen_rand(NUM_RAND_FRAGMENTS, PSF_WIDTH_PIXELS, MEAN_FRAGMENT_LENGTH, FRAGMENT_LENGTH_STD, ...
    ADDED_NOISE_MEAN,ADDED_NOISE_STD,FRAGMENT_STRETCH_STD,IS_CIRCULAR,MINIMUM_LENGTH,TOTAL_RAND_LENGTH,...
    temp_bar);
    % Generate synthetic data for validation
    % version of chrom_assembly_script_vs_ref_synthetic


    % Args:
    %
    % NUM_RAND_FRAGMENTS = 20;% number of random fragments
    % PSF_WIDTH_PIXELS = 4.7; % psf
    % MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
    % FRAGMENT_LENGTH_STD = 100; %std of length of fragment
    % ADDED_NOISE_MEAN = 0.4; % additive noise mean
    % ADDED_NOISE_STD = 0.1; % additive noise std
    % FRAGMENT_STRETCH_STD = 0;%0.02; % fragment length-rescale factor std 0.02
    % IS_CIRCULAR = 1; % whether circular barcode
    %    
    % Returns:
    %   barcodeGen,synthStr,synthStr2,theoryStruct,refBarcode, origPos, origFlip,origStr  
    %   

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

%% STEP1: gen synthetic data

%% gen synthetic
import Nullmodel.gen_random_fragments;
% parameters
% TOTAL_RAND_LENGTH = 10000; % total length

if nargin < 2
    NUM_RAND_FRAGMENTS = 1;
end
if nargin < 2
    PSF_WIDTH_PIXELS = 4.7; % psf
    % NUM_RAND_FRAGMENTS = 20; % number of random fragments
    MEAN_FRAGMENT_LENGTH = 850; % mean length of fragment
    FRAGMENT_LENGTH_STD = 100; %std of length of fragment
    ADDED_NOISE_MEAN = 0.4; % additive noise mean
    ADDED_NOISE_STD = 0.1; % additive noise std
    FRAGMENT_STRETCH_STD = 0;%0.02; % fragment length-rescale factor std 0.02
    IS_CIRCULAR = 1; % whether circular barcode
    MINIMUM_LENGTH = 150;
    TOTAL_RAND_LENGTH = 10000;
end

% if theory, generate
%     
% % psffac = 1;
% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta'};
% %,'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta','sequence1.fasta','sequence2.fasta','sequence3.fasta'};
% nmPerPx = 110;
% nmbp = 0.225;
% psffac = 1/300; % bp level
% % generate theoretical
% nmbp = nmbp*psffac;
% nmPerPx = nmPerPx*psffac;
% import Thry.gen_theoretical;
% [theoryStruct] = gen_theoretical(fastaFile,nmbp,0,nmPerPx);
%     
% temp_bar = importdata(theoryStruct{1}.filename);
    
if nargin >=11
% generate barcodes /also add pcc?
[barcodeGen, refBarcode, origPos, origFlip,origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH,...
    FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD,IS_CIRCULAR, temp_bar);
else
[barcodeGen, refBarcode, origPos, origFlip,origStr] = gen_random_fragments(TOTAL_RAND_LENGTH, PSF_WIDTH_PIXELS, NUM_RAND_FRAGMENTS, MEAN_FRAGMENT_LENGTH,...
    FRAGMENT_LENGTH_STD, ADDED_NOISE_MEAN, ADDED_NOISE_STD, FRAGMENT_STRETCH_STD,IS_CIRCULAR); % random theory
end
% 
% ix=4;
% quick_plot(ix,refBarcode,barcodeGen,origPos(ix), origFlip(ix)+1,1./origStr(ix))

% plot_quick(6,refBarcode,barcodeGen,origPos, origFlip,1./origStr)

% only keep barcodes longer than minimum length.
% minLen = 300; %150kb? or less? depends on application
% barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
% barcodeGen = barcodeGen(barLens>minLen);
% lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);


% plot the fragment positions. Requires hca function
synthStr = cell(1,length(barcodeGen));
for i=1:length(barcodeGen)
    synthStr{i}.idx = 1;
    synthStr{i}.pos = origPos(i);
    synthStr{i}.maxcoef = 0.99; % calculate max coef against thry
    synthStr{i}.lengthMatch = round(length(barcodeGen{i}.rawBarcode)/origStr(i));
    synthStr{i}.lengthUnrescaled = length(barcodeGen{i}.rawBarcode);
    synthStr{i}.bestBarStretch = 1./origStr(i);
    synthStr{i}.or = double(origFlip(i))+1;
    synthStr{i}.rf = origStr(i);
end

% if circ, check if loop around
synthStrLin = synthStr;
% for i=1:length(barcodeGen)
%     if synthStr{i}.pos+synthStr{i}.lengthMatch-1>TOTAL_RAND_LENGTH
%         synthStrLin{i}.pos = synthStrLin{i}.pos-TOTAL_RAND_LENGTH;
%     end
% end
% expected overlap positions
%
synthStr2 = [];
for i=1:length(barcodeGen)
    for j=1:length(barcodeGen)
        if i~=j
            synStemp = synthStrLin;
            % make sure that synStemp{i}.pos-synStemp{j}.pos close to each
            % other
            
            % if distance greater than total length, move the barcode which
            % is more to the right, to the left
            if synStemp{i}.pos>synStemp{j}.pos
                idx = i;
            else
                idx = j;
            end
            if abs(synStemp{i}.pos-synStemp{j}.pos)>abs(abs(synStemp{i}.pos-synStemp{j}.pos)-TOTAL_RAND_LENGTH)
                synStemp{idx}.pos = synStemp{idx}.pos-TOTAL_RAND_LENGTH;
            end

                 
            rF1 = (synthStr{j}.rf); % rescale A to B length
            rF2 = 1;%1/origStr(j);
            
% 

            % length with proper length re-scaling where B is not
            % re-scaled and A re-scaled to B
            lpA =  round(synStemp{i}.lengthMatch*rF1);
            lpB =  round(synthStr{j}.lengthUnrescaled);
         
            synthStr2(i,j).posB = 1;
            synthStr2(i,j).posA = (synStemp{i}.pos-synStemp{j}.pos);

            if synthStr{j}.or==2 % should have j always with or=0 and then switch, or allow it to be flipped?
                pT  = (synStemp{j}.pos+synStemp{j}.lengthMatch-synStemp{i}.pos-synStemp{i}.lengthMatch)+1; % circular?
                synStemp{i}.or =  -synStemp{i}.or+3;
            else
                pT = synthStr2(i,j).posA+1; % 0 should mean identical position
            end
%                 synthStr2(i,j).posA = round((pT)*origStr(j)); % start position is adjusted by origStr(j)
%             else
%%
                % depends if pos or neg.
            if pT > 0
                synthStr2(i,j).posA =round(pT*origStr(j)); % start position is adjusted by origStr(j) to get correct scaling
            else
                synthStr2(i,j).posA =round(pT*origStr(j)); % correct by 1 if before was from 0 
            end
            
          % consider the jth barcode as always at 1 (posB)

            
            pA =  synthStr2(i,j).posA ;
            pB =  synthStr2(i,j).posB ;
            
            st = max(pA,pB)-pB+1;%-( pA-pB+1); % start along un-rescaled reference
            stop = min(pA+lpA-1,pB+lpB-1)-pB+1;
            
      
            re = barcodeGen{i}.rawBarcode;
            re = interp1(re, linspace(1, length(re), round(length(re)*rF1*synthStr{i}.bestBarStretch)));
%             reb = barcodeGen{i}.rawBitmask;
%             reb = interp1(reb, linspace(1, length(reb), round(length(reb)*rF1)));
            rt = barcodeGen{j}.rawBarcode;
            rt = interp1(rt, linspace(1, length(rt), round(length(rt)*rF2)));

            a = [ re re];
            b =[ rt rt];
            if synStemp{i}.or==2
                a = fliplr(a);
            end
%             if synthStr{j}.or==1 % should have j always with or=0 and then switch, or allow it to be flipped?
%                 b = fliplr(b);
%             end
            aFul = a(st-synthStr2(i,j).posA+1:stop-synthStr2(i,j).posA+1);
            bFul = b(st:stop);
            
            % position of overlap
            synthStr2(i,j).pA = st-synthStr2(i,j).posA+1 ;
            synthStr2(i,j).pB = st ;
% 
%             figure,plot(zscore(aFul))
%             hold on
%             plot(zscore(bFul))
%             legend({'i','j'})
            %%
            
%             figure,plot(-synthStr2(i,j).pA:-synthStr2(i,j).pA+length(re)-1,re)
%             hold on
%             plot(-synthStr2(i,j).pB:-synthStr2(i,j).pB+length(rt)-1,rt)

% %       
            synthStr2(i,j).fulloverlapPosA = pB-st+1:pB+stop-1;
            synthStr2(i,j).fulloverlapPosARescaled = pA-st+1:pA+stop-1;
            synthStr2(i,j).fullscore = zscore(aFul(10:end-9),1)*zscore(bFul(10:end-9),1)'/length(aFul(10:end-9));

            synthStr2(i,j).overlaplen = length(aFul);
            if synthStr2(i,j).overlaplen >= MINIMUM_LENGTH
                synthStr2(i,j).score =   synthStr2(i,j).fullscore ;
            else
               synthStr2(i,j).score = nan;
            end
            synthStr2(i,j).bestBarStretch = rF1*synthStr{i}.bestBarStretch;
            synthStr2(i,j).or = synStemp{i}.or;
            
            synthStr2(i,j).lenB = length(barcodeGen{j}.rawBarcode);
            synthStr2(i,j).lenA = length(interp1(barcodeGen{i}.rawBarcode, linspace(1, length(barcodeGen{i}.rawBarcode), round(length(barcodeGen{i}.rawBarcode)*synthStr2(i,j).bestBarStretch))));
%             synthStr2(i,j).fullscore = zscore(aFul,1)*zscore(bFul,1)'/length(aFul);
%             overlapStruct(k,iy).overlaplen = length(aFul);
%             overlapPos = 
        else
            synthStr2(i,j).fullscore =nan;
            synthStr2(i,j).score =nan;

        end   
    end
end
% nd

[a,msg] = mkdir(strcat('output',timestamp));
i=1;theoryStruct=[];
theoryStruct{i}.filename = fullfile(strcat('output',timestamp),'theory_barcode.txt');
fileID = fopen(theoryStruct{i}.filename,'w');
fprintf(fileID,strcat([' %2.' num2str(14) 'f ']), refBarcode);
fclose(fileID);
theoryStruct{i}.meanBpExt_nm = 0.3;
theoryStruct{i}.psfSigmaWidth_nm = 600;
theoryStruct{i}.length = length(refBarcode);
theoryStruct{i}.isLinearTF = 0;
theoryStruct{i}.name = 'Synthetic theory';

sets.comparisonMethod = 'mass_pcc';
sets.filter = 0;
sets.filterSettings.filter = 0;
% nmbp = barcodeGen{1}.nmbp;
nmpx = 110; % 208?
psf = 300; %300 nm

end