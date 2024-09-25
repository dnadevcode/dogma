
%
% Reference based validation

[outConsensus, coverage, pval] = gen_reference_based_assembly(barcodeGen,synthStr{idxRun},theoryStruct{idxRun},timestamp);

bar = nanmean(outConsensus);

cGen =[];
 cGen{1}.rawBarcode = bar;

 cGen{1}.rawBitmask = ~isnan(bar);

w = 500;% minimum overlap length
sF = 0.9:0.025:1.1;

bT{1}.rawBarcode = theoryStruct{idxRun}{1}.rawBarcode;
bT{1}.rawBitmask = ones(1,length(bT{1}.rawBarcode));
bT{1}.rawBitmask = logical(bT{1}.rawBitmask);

wBoostrapping = 1000;
numWindows = floor(length( cGen{1}.rawBarcode)/wBoostrapping);
R = reshape(1:floor(length(cGen{1}.rawBarcode)/numWindows)*numWindows,floor(length(cGen{1}.rawBarcode)/numWindows),[]);
     

% todo: how much does this improve local scores? We can draw randomly (the
% same as for synthetic barcodes), and then compare the distributions.
% Should be improved for drawing from a averaged one
f = figure;
tiledlayout(10,2,'TileSpacing','tight');

for k=1:10;
    cos{1}.rawBarcode = cGen{1}.rawBarcode(R(:,k));
    cos{1}.rawBitmask = logical(cGen{1}.rawBitmask(R(:,k)));
    
    %  add theory as last barcode
    consensusGen =[cos bT];
    
    % calculate MP scores 
    [oScons{1}] = calc_overlap_mp(consensusGen,sF, w,timestamp);
     
    % [oScons] = calc_scores({cGen},w,sF);
    
    
           
    %% Plot specific pair
    
    import Core.plot_match_simple;
    expidx = 1; % idx of experiment
    thrIdx = 2; % idx of theory
    lenThry = theoryStruct{idxRun}{1}.length;
    import Plot.pair_evaluation_with_ground_truth_plot;
    [f] = pair_evaluation_with_ground_truth_plot(consensusGen, oScons{1},expidx,thrIdx,[],lenThry,f);
end
