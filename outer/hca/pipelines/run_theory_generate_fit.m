function [t,theoryGen] = run_theory_generate_fit(fold,nmbp,yoyoBindingProb,precalc)

% Inline theory generation without passing through hca_barcode_alignment for the fitted YOYO-1 binding probability model,
% One can use this function with batch command
%  c = parcluster;
% batch(c,@run_theory_generate_fit,1,{folderName,nmbp,0},'Pool',29);

t0 = tic;

% addpath(genpath('/home/avesta/albertas/reps/hca'));

timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');

if precalc==1
    a=dir(fullfile(fold,'*.mat'));
    folderName = arrayfun(@(x) fullfile(a(x).folder,['seq', num2str(x),'.mat']),1:length(a),'un',false);
else
    a=dir(fullfile(fold,'*.fasta'));
    folderName = arrayfun(@(x) fullfile(a(x).folder,a(x).name),1:length(a),'un',false);
end


[hcatheory.sets,hcatheory.names] = Core.Default.read_default_sets('hcasets.txt');

hcatheorySets = hcatheory.sets.default;

hcatheorySets.folder = folderName;

hcatheorySets.computeBitmask = 0;
if nargin>=2
    hcatheorySets.meanBpExtNm = nmbp;
end

disp(['N = ',num2str(length(  hcatheorySets.folder )), ' sequences to run'])

theoryGen = struct();

% save sets
theoryGen.sets = hcatheorySets;

% sigma = 0.68; % scaling of AT-GCs
gcSF  = 1;
% kY = 10; % binding constants, hard-coded for now
% kN = 30;
psf = hcatheorySets.psf; % seemed to give best results
ligandLength = hcatheorySets.ligandLength;
% cN = hcatheorySets.concN;
% cY = hcatheorySets.concY;
nmpx = hcatheorySets.pixelWidthNm;
isC = ~hcatheorySets.isLinearTF;
pxSize = hcatheorySets.pixelWidthNm/hcatheorySets.meanBpExtNm;

folds = hcatheorySets.folder;


theoryBarcodes = cell(1,length(folds));
theoryBitmasks = cell(1,length(folds));
theoryNames = cell(1,length(folds));
theoryIdx = cell(1,length(folds));

% parlist = [gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY,cN,kY,ligandLength];
% parlistcell = num2cell(parlist) ;

% C = parallel.pool.Constant(parlistcell);

%todo: add a progressbar
import CBT.SimpleTwoState.gen_simple_theory_px_fit;
        import CBT.SimpleTwoState.px_cut_pos;



% N = length(folds);
% parfor_progress(N); % Initialize 

parfor idx=1:length(folds)
%     idx
    data = struct("atsum",[],'name',[],'idsElt',[],'pxcut',[]);

    if precalc
       data = load(folds{idx});
    else
    % slow part: reading sequence and indexing. Could be pre-calculated
        fasta = fastaread(folds{idx});
        data.ntSeq = nt2int(fasta.Sequence,'Unknown',4);
        data.ntSeq(data.ntSeq>4) = 4; % simplify by converting all special letters to thymidine. Works fine if very few NN
    
        sz = ones(1,ligandLength)*4;
        I = arrayfun(@(x) data.ntSeq(x+(1:length(data.ntSeq)-ligandLength+1)),0:ligandLength-1,'un',false);
        data.idsElt = sub2ind(sz, I{:} );
    
        data.atsum = cumsum((data.ntSeq == 1)  | (data.ntSeq == 4) );
        data.name = fasta.Header;

        [data.pxcut.pxCutLeft, data.pxcut.pxCutRight, data.pxcut.px] = px_cut_pos( data.atsum, gcSF, pxSize);
 
    end
    %

    % fast part
    [theorySeq] = gen_simple_theory_px_fit(data.pxcut,gcSF,pxSize,nmpx,isC,[],[],psf,[], [],[],ligandLength,yoyoBindingProb,data.idsElt);


%     [theorySeq] = gen_simple_theory_px(numWsCumSum,gcSF,pxSize,nmpx,isC,sigma,kN,psf,cY, cN,kY);


    theoryBarcodes{idx} = theorySeq;
    theoryBitmasks{idx} = [];
    
    theoryNames{idx} = data.name;
    theoryIdx{idx} = idx;
%     if mod(idx,200) == 0
%         parfor_progress; % Count 
%     end
%     toc
end
%    parfor_progress(0); % Clean up

disp('Finished calculating theories')


    
    theoryGen.theoryBarcodes = theoryBarcodes;
    theoryGen.theoryBitmasks = theoryBitmasks;
    theoryGen.theoryNames = theoryNames;
    theoryGen.theoryIdx = theoryIdx;
%     theoryGen.bpNm = bpNm;


    

    resultsDir = fullfile(fileparts(folds{1}),'theoryOutput'); % we always save theories output in the same folder as provided data
    
    
    % make theoryData folder
    [~,~] = mkdir(resultsDir);


    matFilename = strcat(['theoryGen_', num2str(hcatheorySets.meanBpExtNm) '_' num2str(hcatheorySets.pixelWidthNm) '_' num2str(psf) '_' num2str(~isC) '_' sprintf('%s_%s', timestamp) 'session.mat']);
    matFilepath = fullfile(resultsDir, matFilename);
    

%     if nargout < 1
%         assignin('base','theoryGen', theoryGen)
%         disp('Assigned theoryGen to workspace');
        
        
        save(matFilepath, 'theoryGen','-v7.3');
        fprintf('Saved theory mat filename ''%s'' to ''%s''\n', matFilename, matFilepath);
% %     end


% import Core.run_hca_theory;
% run_hca_theory(hcatheorySets);

t = toc(t0);

end



