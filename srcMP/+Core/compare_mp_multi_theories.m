function [compStrOut,mp1individualOut] = compare_mp_multi_theories(barcodeGen,theoryStruct,sF,MIN_OVERLAP_PIXELS_BL,numWorkers)
% like compare_mp_all but for multi theories. Using SCAMP if available,
% alternatively local MP function

if nargin < 5
    numWorkers = 30;
end
% MIN_OVERLAP_PIXELS_BL = 150;
MIN_OVERLAP_PIXELS = MIN_OVERLAP_PIXELS_BL;
% minLen = 500; %150kb? or less? depends on application // if GUI, user selects thi
timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
out=strcat('output',timestamp);
mkdir(out);


% barLens = cellfun(@(x) sum(x.rawBitmask),barcodeGen);
% barcodeGen = barcodeGen(barLens>MIN_OVERLAP_PIXELS_BL);
lengths = cellfun(@(x) sum(x.rawBitmask),barcodeGen);

barStruct = cell2struct([cellfun(@(x) x.rawBarcode,barcodeGen,'un',false);...
    cellfun(@(x) x.rawBitmask,barcodeGen,'un',false)]',{'rawBarcode','rawBitmask'},2);

% 
% if ispc % for PC get from the initialization
%     sets.SCAMP_LINE = 'C:\Users\Lenovo\git\SCAMP' ;
%     SCAMP_LINE = strcat([sets.SCAMP_LINE '\build\Release\SCAMP.exe']);  %'.\..\SCAMP\build\Release\SCAMP.exe'; % windows
% else %
% %     sets.SCAMP_LINE = '/home/albyback/postdocData/test_transloc/SCAMP/';
%     sets.SCAMP_LINE = '~/SCAMP/';
% 
%     SCAMP_LINE = strcat([sets.SCAMP_LINE '/build/SCAMP']); %'~/git/SCAMP/build/SCAMP'; % how to call scamp from bash
% end
SCAMP_LINE = 'SCAMP';


foldSynth = 'barcoli';

% [namesBar, stridx] = Core.save_bars_rescaled_txt(barStruct,sF,foldSynth);
[namesBar, stridx, baridx] = Core.save_bars_rescaled_txt_stack(barStruct,sF,foldSynth);
stridx = stridx(1:end-MIN_OVERLAP_PIXELS+1);
baridx = baridx(1:end-MIN_OVERLAP_PIXELS+1);

import Core.calc_overlap;
sFList = [1:length(sF) -(1:length(sF))];


for ix=1:length(theoryStruct)
% names2{1} = ;// could also stack the two theories (esp only if we're
% interested in best match)
    try
        names2 = arrayfun(@(x) theoryStruct{ix}.filename, 1:length(namesBar),'un',false);
    catch
        names2 = arrayfun(@(x) theoryStruct(ix).filename, 1:length(namesBar),'un',false);
    end

% names2 = arrayfun(@(x) theoryStruct{x}.filename, 1:length(theoryStruct),'un',false);

%     tic
    [mpI1,mp1,maxMP] = calc_overlap([],timestamp,SCAMP_LINE, MIN_OVERLAP_PIXELS, numWorkers,namesBar,names2);
%     toc

% 
%     % % %% Check the re-scaling factors along the barcode. 
%     strFac = cell(1,length(barStruct));
%     pos = cell(1,length(barStruct));
%     bestSTF = zeros(1,length(barStruct));
    compStr = cell(1,length(barStruct));
%     mp1individual = cell(1,length(barStruct));
%     mpI1individual = cell(1,length(barStruct));
%     stridxindividual = cell(1,length(barStruct));
%     for i=1:length(barStruct)
%         mp1individual{i} = mp1{1}(baridx==i);
%         mpI1individual{i} = mpI1{1}(baridx==i);
%         stridxindividual{i} = stridx(baridx==i);
%     
%         locs = [];
%         meanMP=[];
%         maxMp = [];
%         for j=1:length(sFList) % length re-scaling based on average of a few highest values.
%             mpS = mp1individual{i}(stridxindividual{i}==sFList(j));
%             [mpS,locs{j}] = sort(mpS,'descend','MissingPlacement','last');
%             meanMP(j) = nanmean(mpS(1:min(end,20))); %1 would be standard
%             maxMp(j) = mpS(1);
%             locs{j} = locs{j}(~isnan(mpS));
%         end
% %         compStr{i}.meanMP = meanMP;
%         [aval,bstr] = max(meanMP); % find best re-scale factor
%     
%         strFac{i} = sFList(bstr); % +1??
%     
%         pos{i} = locs{bstr}(1);
%     
%         mpIs = mpI1individual{i}(stridxindividual{i}==sFList(bstr));
%     
%     
%         bestSTF(i)=  sF(abs(strFac{i}(1)));
%         compStr{i}.pB = mpIs(locs{bstr}(1))'+1;
%         compStr{i}.pos = mpIs(locs{bstr}(1))'+1;
%     
%     
%         compStr{i}.pA =  pos{i};% position on strFac - from this can extract overlap pos 
%         compStr{i}.or = sign( strFac{i}); % should be able to plot exact match..
%         compStr{i}.idx = ix;
%         compStr{i}.bestBarStretch =   bestSTF(i);
%         compStr{i}.lengthMatch =  lengths(i)*bestSTF(i);
%         compStr{i}.maxcoef = maxMp(bstr)';
% %         compStr{i}.allposA = locs{bstr};
% %         compStr{i}.allposB = mpIs(locs{bstr})+1;
%     
%     end

    compStrOut{ix} = compStr;
    if nargout > 1
        mp1individualOut{ix} = mp1individual
    end
end

end

