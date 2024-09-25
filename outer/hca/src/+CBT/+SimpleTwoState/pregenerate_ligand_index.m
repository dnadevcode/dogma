function [] = pregenerate_ligand_index(fold,fold2,ligandLength,gcSF, pxSize)
% Pre-generates index for ligand for faster theory generation

a = dir(fullfile(fold,'*.fasta'));
folderName = arrayfun(@(x) fullfile(a(x).folder,a(x).name),1:length(a),'un',false);
    
import CBT.SimpleTwoState.px_cut_pos;



parfor idx=1:length(folderName)
    % idx

    % slow part: reading sequence and indexing. Could be pre-calculated
    

    fasta = fastaread(folderName{idx});
    ntSeq = nt2int(fasta.Sequence,'Unknown',4);
    ntSeq(ntSeq>4) = 4; % simplify by converting all special letters to thymidine. Works fine if very few NN

    sz = ones(1,ligandLength)*4;
    I = arrayfun(@(x) ntSeq(x+(1:length(ntSeq)-ligandLength+1)),0:ligandLength-1,'un',false);
    idsElt = sub2ind(sz, I{:} );

    atsum = cumsum((ntSeq == 1)  | (ntSeq == 4) );
    

% OLD CODE  - slow  
%     tic
%     idsElt = zeros(1,length(ntSeq)-ligandLength+1);
%     for i=1:length(ntSeq)-ligandLength+1
%         cellInd = num2cell(ntSeq(i:i+ligandLength-1));
%         idsElt(i) = sub2ind(sz, cellInd{:} );
%     end
%     toc
%     isequal(idsElt2,idsElt)

    pxcut = struct("pxCutLeft",[],"pxCutRight",[],"px",[]);
    [pxcut.pxCutLeft, pxcut.pxCutRight, pxcut.px] = px_cut_pos( atsum, gcSF, pxSize);
       

    name = fasta.Header;

    s = struct("atsum",atsum,'name',name,'idsElt',idsElt,'pxcut', pxcut);
%     save(sprintf("output_%d.mat",idx),"-fromstruct",s);

    save(fullfile(fold2,['seq', num2str(idx),'.mat']),"-fromstruct",s);

% 
%     tic 
%     load(fullfile('/export/scratch/albertas/download_dump/single2/',['seq', num2str(1),'.mat']))
%     toc
end

end

