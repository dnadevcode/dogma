lenSeq = 5.09*10^6;
nmPerPx = 110; % there's some dependence on these
nmbp = 0.2; % nm/bp (from lambdas
psffac = 1; % scaling factor for psf
numWorkers = 4; % num workers ( for parpool)
minLen = 100:50:2500;

NN = 10;
mpMax = cell(1,NN); 
for i=1:NN
    [mp,mpMax{i}] = bargrouping_minimum_length([]); % maybe need to define length? (to be the same)
end

mpMax = cell2mat(mpMax')
% run for the whole database?
% fastaFile = {"/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/DA32087.fasta"};
fastaFile = {"/proj/snic2022-5-384/users/x_albdv/data/bargrouping/ecoli/FASTAS/018_final_polish.fasta"};
[MP,mpMax,theoryStructRev,MPI,theoryStruct] = bargrouping_minimum_length(fastaFile,nmPerPx,nmbp,psffac,numWorkers,minLen);


% figure, errorbar(minLen,mean(cell2mat(mpMax')),std(cell2mat(mpMax')))
% hold on
figure
plot(minLen,mpMax)
xlabel('Overlap length')
ylabel('Max overlap PCC')
legend({'E-coli'})

% legend({'Simulated','E-coli'})

%% Plot individual

ii=1;
scores = MP{ii}{1}(1:end/2);
[a,b] = findpeaks(scores,'MinPeakDistance',200,'SortStr','descend');

iy = 1;
osTheory = [];
for k=b'
    osTheory(k,iy).pA = mod(MPI{ii}{1}(k),theoryStruct{1}.length)+1;
    % if greater than barcode length, then it's the reverse
    osTheory(k,iy).pB = k;
    osTheory(k,iy).h = minLen(ii);
    osTheory(k,iy).or =  sign(theoryStruct{1}.length -MPI{ii}{1}(k));
    osTheory(k,iy).bestBarStretch = 1; % no rescaling;
    osTheory(k,iy).score = scores(k); % no rescaling;
    osTheory(k,iy).fullscore = nan; % no rescaling;

end
barStruct = [];
barStruct(1).rawBarcode = data(2:end/2);
barStruct(1).rawBitmask = logical(ones(1,length(barStruct(1).rawBarcode )));
% 
% [sortSc,sortId] = sort(scores,'descend','MissingPlacement','last');
%% Instead use the other function
ix=b(2);

barStruct(ix) = barStruct(1);
import Core.plot_match_simple;
% [f] = plot_match_simple(barStruct, oS,curSink,curSource);
[f] = plot_match_simple(barStruct, osTheory,ix,1);
