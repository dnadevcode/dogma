%% This script runs bargrouping_w_theory without hpc, i.e. reference based chromosome assembly. HPC commands are commented out (%)

%% this would run on kebnekaise.
% c = parcluster('kebnekaise');
% c.AdditionalProperties.AccountName = 'snic2021-5-132';
% c.AdditionalProperties.WallTime = '12:00:00';
% batchPoolSize = 119;

%
cmap = colormap;
numc = length(cmap);
close(gcf);

%% General parameters

MIN_OVERLAP_PIXELS = 100; % minimum number of overlap between two barcodes
STRETCH_MAX = 0.1; % maximum stretch
STRETCH_STEP = 0.005; % stretch step
MIN_BARCODE_SCORE_THRESHOLD = 3;
MIN_BARGROUP_SCORE_THRESHOLD = 7; % threshold for belonging to bargroup 
MIN_REFGROUP_SCORE_THRESHOLD = 0; % threshold for belonging to refgroup (so 0 means everything is mapped somewhere)
P_COMBINE_METHOD = 'stouffer';
PRINT_TO_WINDOW = 1;

NMPX = 130; % experiment parameters, nm/px camera parameter
NMPSF = 300; % psf in nanometers
NMBP = 0.28; % nm/bp extension, typically if many different experiments, average nm/bp best

numUncert = ceil(5*NMPSF/NMPX); % uncertainty in pixels on left and right edge of barcode

%% INPUT DATA

mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/Yeast 280721'; % main path with kymo's. Change this to yours
savepath = strcat('output_',datestr(datetime)); % create output directory for saving results
mkdir(savepath);
kymoPath = {dir([mainPath '/Raw kymos/*.tif'])}; % extract kymo's from this folder
kymoParams = [];

% link to DNA SEQUENCE
seqPath = fullfile(mainPath, 'sequence.fasta');


% TO BEGIN, maybe also run the following two lines, to see how much time it
% takes  with few barcodes (i.e. 50), and low threshold (i.e. MIN_BARGROUP_SCORE_THRESHOLD=4)
% kymoPath{1} = kymoPath{1}(1:min(end,50));
% MIN_BARGROUP_SCORE_THRESHOLD = 4; % threshold for belonging to bargroup 


% a number of commented out test folders.

% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/DA32087, 180815';
% kymoPath = {dir([mainPath '/RawKymographs 180815 - 130nm - 0.273 nm_bp/*.tif'])};
% kymoParams = [];

% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/DA65787, 190710_';
% kymoPath = {dir([mainPath '/Raw Kymographs 190710 - 130nm - 0.274nm_bp/*.tif']);
%     dir([mainPath '/Raw Kymographs 190710 - 208nm - 0.274nm_bp/*.tif'])};
% kymoParams = struct('nmpx', [130 208], 'nmbp', [0.274 0.274]);

% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/DA69401, 200708';
% kymoPath = {dir([mainPath '/1,6 x Fragments Kymos 0.239/*.tif']);
%     dir([mainPath '/1 x Fragments Kymos 0.239/*.tif'])};
% kymoParams = struct('nmpx', [130 208], 'nmbp', [0.239 0.239]);

% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/DA69402_200630';
% kymoPath = {dir([mainPath '/1,6 x Fragments Kymos 0.281/*.tif']);
%     dir([mainPath '/1 x Fragments Kymos 0.281/*.tif'])};
% kymoParams = struct('nmpx', [130 208], 'nmbp', [0.281 0.281]);
% 
% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/DA69405, 200630_';
% kymoPath = {dir([mainPath '/1,6 x Fragments Kymos 0.287/*.tif']);
%     dir([mainPath '/1 x Fragments Kymos 0.287/*.tif'])};
% kymoParams = struct('nmpx', [130 208], 'nmbp', [0.287 0.287]);

% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/DA69406, 200630';
% kymoPath = {dir([mainPath '/1,6 x Fragment Kymos 0,264/*.tif']);
%     dir([mainPath '/1 x Fragment Kymos 0,264/*.tif'])};
% kymoParams = struct('nmpx', [130 208], 'nmbp', [0.264 0.264]);

% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/EF512 (Sample 4)';
% kymoPath = {dir([mainPath '/RawKymographs 200527 - 130nm - 0.282nm_bp/*.tif']);
%     dir([mainPath '/RawKymographs 200527 - 208nm -0.282nm_bp/*.tif']);
%     dir([mainPath '/RawKymographs 20200324 - 160nm - 0.232nm_bp/*.tif'])};
% kymoParams = struct('nmpx', [130 208 160], 'nmbp', [0.282 .282 0.232]);

% mainPath = '/proj/nobackup/snic2021-5-132/data/Chromosome/EF575 (Sample 9)';
% kymoPath = {dir([mainPath '/RawKymographs 130nm - 200529 - 0.272nm_bp/*.tif']);
%     dir([mainPath '/RawKymographs 160nm - 20200422 - 0.3nm_bp_/*.tif']);
%     dir([mainPath '/RawKymographs 160nm - 20200422 - 0.268nm_bp_/*.tif']);
%     dir([mainPath '/RawKymographs 208nm - 200529 - 0.272nm_bp/*.tif'])};
% kymoParams = struct('nmpx', [130 160 160 208], 'nmbp', [0.272 .3 0.268 0.272]);


%% load and align barcodes

import TempStuff.gen_barcodes_from_kymos % this extracts and aligns barcodes. TODO: change this to use HCA to maybe run less time-frames or simpler alignment method
barcodeGen = gen_barcodes_from_kymos(kymoPath, NMPSF/NMPX);
import TempStuff.rescale_barcodes
if not(isempty(kymoParams))
  barcodeGen2 = rescale_barcodes(barcodeGen, kymoParams, NMPX, NMBP);
else
  barcodeGen2 = barcodeGen;
end
barcodeGen3 = vertcat(barcodeGen2{:});

%%

barcodes = {barcodeGen3.rawBarcode}';
bitmasks = cellfun(@(x) [false(1, numUncert) true(1, size(x, 2) - 2 * numUncert) false(1, numUncert)], barcodes, 'un', 0);
bitmasks = cellfun(@(x, y) y & abs(x - mean(x(y))) <= 2.5 * std(x(y)), barcodes, bitmasks, 'un', 0);
% bitmasks = cellfun(@(x) true(size(x)), barcodes, 'un', 0);

%% generates theory sequence

fasta = fastaread(seqPath);
ntSeq = nt2int(fasta.Sequence);
import Theory.compute_theory_barcode
refBarcode = compute_theory_barcode(ntSeq, NMPX/NMBP, NMPSF/NMBP, true);
refBitmask = movsum(ntSeq > 4, ceil(NMPX/NMBP)) < 50;
refBitmask = refBitmask(round(linspace(1, length(ntSeq), length(refBarcode))));

import Core.Bargroup
refId = length(barcodes)+1;
refBargroup = Bargroup(refId, 1, 0, 1, length(refBarcode));

%% Run all to all comparison for all re-scaling factors


[scoreRegistry,scoreRegistryIndex] = Core.compute_all_barcode_pair_scores([barcodes; refBarcode], [bitmasks; refBitmask], MIN_BARCODE_SCORE_THRESHOLD, STRETCH_MAX, STRETCH_STEP, MIN_OVERLAP_PIXELS, P_COMBINE_METHOD, PRINT_TO_WINDOW);

%% HPC version
% j = c.batch(@Core.compute_all_barcode_pair_scores, 2, {[barcodes; refBarcode], [bitmasks; refBitmask], MIN_BARCODE_SCORE_THRESHOLD, STRETCH_MAX, STRETCH_STEP, MIN_OVERLAP_PIXELS, P_COMBINE_METHOD, 0}, 'Pool', batchPoolSize);
% j.wait;
% out = j.fetchOutputs();
% scoreRegistry = out{1};
% scoreRegistryIndex = out{2};

%% Create bargroups
[outputTree,~, mergeScores] = Core.bargroup_hc(barcodes, scoreRegistry, scoreRegistryIndex, 0, STRETCH_MAX, STRETCH_STEP, P_COMBINE_METHOD,MIN_OVERLAP_PIXELS, PRINT_TO_WINDOW);

%% hpc version
% j = c.batch(@Core.bargroup_hc, 3, {barcodes, scoreRegistry, scoreRegistryIndex, 0, STRETCH_MAX, STRETCH_STEP, P_COMBINE_METHOD, MIN_OVERLAP_PIXELS, PRINT_TO_WINDOW}, 'Pool', batchPoolSize, 'CaptureDiary', true);
% j.wait;
% out = j.fetchOutputs();
% outputTree = out{1};
% mergeScores = out{3};

%% plot
figure; plot(mergeScores); %hold on; plot(cellfun(@(x) std(cellfun(@(y) length(y.ids), x)), outputTree))
savefig(fig, fullfile(savepath, 'mergeScores'));

%% CREATE bagroups

bargroups = outputTree{find(mergeScores < MIN_BARGROUP_SCORE_THRESHOLD, 1, 'first') - 1};%end};%

%% --- Plot results ---
% here we can change 
MIN_REFGROUP_SCORE_THRESHOLD
% in order to get different output figures

fprintf('\nFinal results:\n')

mergedBargroup = refBargroup;
bargroupScore = zeros(1, length(bargroups));
rootOffset = zeros(1, length(bargroups));
isFlipped = zeros(1, length(bargroups));
stretchAmount = zeros(1, length(bargroups));
cmapIdx = ceil(randperm(length(bargroups))*numc/length(bargroups));
bgColors = zeros(length(barcodes), 3);
memcount = cellfun(@(x) length(x.ids), bargroups);

for i = 1:length(bargroups)
%     if length(bargroups{i}.ids) > 1
        [bargroupScore(i), rootOffset(i), isFlipped(i), stretchAmount(i)] = refBargroup.compare(bargroups{i}, scoreRegistry, scoreRegistryIndex, STRETCH_MAX, STRETCH_STEP, P_COMBINE_METHOD, false);
        if bargroupScore(i) >= MIN_REFGROUP_SCORE_THRESHOLD
            mergedBargroup = mergedBargroup.merge(bargroups{i}, rootOffset(i), isFlipped(i), stretchAmount(i), bargroupScore(i));
            bgColors(bargroups{i}.ids, :) = repmat(cmap(cmapIdx(i), :), length(bargroups{i}.ids), 1);
    %     elseif length(bargroups{i}.ids) > 2
    %         disp('a')
    %         bargroups{i}.visualize([barcodes; refBarcode], false, true);
        end
%     end
end
fig = mergedBargroup.visualize([barcodes; refBarcode], false, true, setxor(mergedBargroup.ids, refId), refId, 'Reference', [bitmasks; refBitmask], bgColors);
finalScore = refBargroup.compare(mergedBargroup, scoreRegistry, scoreRegistryIndex, STRETCH_MAX, STRETCH_STEP, P_COMBINE_METHOD, true);
fprintf('Barcodes:%.0f\tBargroups:%.0f\tScore:%.2f\n', length(mergedBargroup.ids) - 1, sum(cellfun(@(x) length(x.ids), bargroups(bargroupScore >= MIN_REFGROUP_SCORE_THRESHOLD)) > 1), finalScore)
savefig(fig, fullfile(savepath, compose("bargrouping_res_T=%g", MIN_REFGROUP_SCORE_THRESHOLD)));
% savefig(fig, fullfile(mainPath, compose("bargrouping_res_only_groups_T=%g", MIN_REFGROUP_SCORE_THRESHOLD)));

%%

mergedGroupIds = mergedBargroup.ids(2:end);
barcodePos = zeros(size(mergedGroupIds));
%
for i = 1:length(mergedGroupIds)
    thisId = mergedGroupIds(i);
    tmpBargroup = Core.Bargroup(thisId, 1, false, 1, length(barcodes{thisId}));
    [~, barcodePos(i)] = refBargroup.compare(tmpBargroup, scoreRegistry, scoreRegistryIndex, STRETCH_MAX, STRETCH_STEP, P_COMBINE_METHOD, true);
end
%
fig = figure;
histogram(abs(mergedBargroup.pos(2:end) - barcodePos(:) - mergedBargroup.pos(1)), 20);
savefig(fig, fullfile(savepath, 'bargrouping_dist'));

%%

save(fullfile(savepath, 'bargrouping_res'), 'barcodes', 'bitmasks', 'scoreRegistry', 'scoreRegistryIndex', 'outputTree', 'mergeScores', 'bargroups', 'refBarcode', 'refBitmask', 'mergedBargroup');