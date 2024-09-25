%%
function [MP,mpMax,theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,psffac,numWorkers,minLen, lenSeq)


% Script to check local similarities of theoretical barcodes in order to
% estimate the minimum overlap length/ 

% In practice has dependece both on minimum length AND the lengths of two
% barcodes, since when we slide one barcode on top of the other, the length ranges from 
% "min length" to "max length". We could limit this "max length" if we want
% to match barcodes with significant overlap, but not "too long" (then it
% appears almos t as a clone of the same position)


% todo: also allow length re-scaling (so do A vs B comparison) 


if isempty(fastaFile)
    simulated = 1;
else
    simulated = 0;
end

% fastaFile =  {'/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/DA32087.fasta','/proj/snic2022-5-384/users/x_albdv/data/E. coli Assembly - data for Lund/018_final_polish.fasta','sequence1.fasta','sequence2.fasta','sequence3.fasta'};

if nargin < 2
%     simulated = 1;
    lenSeq = 5*10^6;
    nmPerPx = 110; % there's some dependence on these
    nmbp = 0.2; % nm/bp (from lambdas
    psffac = 1; % scaling factor for psf
    numWorkers = 4; % num workers ( for parpool)
    % minLen = 200; % min overlap
    minLen = 100:50:2500;
end


if simulated
    [theoryStructRev,theoryStruct,barcodeGen] = prep_thry_for_local_comp_simulated(lenSeq, nmbp, nmPerPx, psffac);

else
    % just DA:
%     fastaFile = {"C:\Users\Lenovo\git\hca\test\features\DA32087.fasta"};
%     fastaFile = {"/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/DA32087.fasta"};
%     fastaFile = {"/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/018_final_polish.fasta"};

    [theoryStructRev,theoryStruct,barcodeGen] = prep_thry_for_local_comp(fastaFile, nmbp, nmPerPx, psffac);

end

% local compare

MP = cell(1,length(minLen));
MPI = cell(1,length(minLen));
for ii=1:length(minLen)
    [MP{ii}, MPI{ii}] = local_compare_thry(minLen(ii),[], nmbp, nmPerPx, psffac, numWorkers,theoryStructRev);
end
mpMax = cellfun(@(x) max(x{1}(1:end/2)),MP);

%% create OS struct
% import Core.mp_res_to_struct;
% [overlapStruct2] = mp_res_to_struct(MP{1},MPI{1},baridx2,stridx,minLen(1),1,barStruct);


% figure,plot(minLen,mpMax)
% xlabel('Overlap length')
% ylabel('Max overlap PCC')


% todo: if we noisify theoryStructRev, we get some error bars for this
% figure,histogram(MP{1}{1}(1:end/2));

% so, if overlap lengths is chosen to be minLen

% [MPthr2, MPIthr2] = local_compare_thry(300);
% MP
% 
% figure
% plot(MPthr1{1})
% hold on
% plot(MPthr2{1})
% xlabel('Position (px)')
% legend({'PCC overlap 150','PCC overlap 300'})
%%
% 
% import Thry.gen_theoretical;
% [theoryStruct,~,barcodeGenT] = gen_theoretical(fastaFile,nmbp,0,nmPerPx);
% 
% sets.comparisonMethod = 'mass_pcc';
% 
% sF = 0.95:0.01:1.05;
% [comparisonStruct] = compare_to_t(bars,theoryStruct,sF,sets)
% 
% figure,histogram(cellfun(@(x) x.bestBarStretch,comparisonStruct))
% 
% mCC=cellfun(@(x) x.maxcoef(1),comparisonStruct);
% mO =cellfun(@(x) x.lengthMatch,comparisonStruct);
% 
% figure,plot(mO,mCC,'x')
% hold on
% plot(minLen,mpMax)
% 
% % we can also say which barcodes are below this line or above
% 
% 
% bpPx = 500;
% import Plot.plot_best_pos;
% fig1 = plot_best_pos([], comparisonStruct, [], [], [],cumsum(cellfun(@(x) x.length,theoryStruct)),1);%1/bpPx*10^6

