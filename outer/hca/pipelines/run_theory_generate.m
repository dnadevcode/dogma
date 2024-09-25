function [t] = run_theory_generate(fold)

% Inline theory generation without passing through hca_barcode_alignment,
% One can use this function with batch command
%  c = parcluster;
% batch(c,@run_theory_generate,1,{folderName},'Pool',2);

a=dir(fullfile(fold,'*.fasta'));
folderName = arrayfun(@(x) fullfile(a(x).folder,a(x).name),1:length(a),'un',false);

t0 = tic;

[hcatheory.sets,hcatheory.names] = Core.Default.read_default_sets('hcasets.txt');

hcatheorySets = hcatheory.sets.default;

hcatheorySets.folder = folderName;

hcatheorySets.computeBitmask = 0;
hcatheorySets.meanBpExtNm = 0.34;

disp(['N = ',num2str(length(  hcatheorySets.folder )), ' sequences to run'])

import Core.run_hca_theory;
run_hca_theory(hcatheorySets);

t = toc-t0;


end


