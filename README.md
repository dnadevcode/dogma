# DOGMA

requires (1) SCAMP CLI for MP based bargrouping. (https://github.com/zpzim/SCAMP)

'''
git clone https://github.com/zpzim/SCAMP
cd SCAMP
git submodule update --init --recursive
mkdir build && cd build
'''

requires (2) most recent version of HCA to compare alignment (https://gitlab.com/dnadev/hca) which is included in subfolder (outer)

requires (3) access to matlab (either desktop installation or matlab online. Tested with MATLAB R2022a)

We provide a live script that generates figures similar to the main figures in the main text
'''
dogma_demo.mlx
'''

paper_figs/final provides scripts used to generate figures in the paper, i.e.
fig1_ex creates Fig1,

final_figure_N (N=3-7) creates figures Fig3-Fig7,

These require additional pre-generated results data that will be put on a public database with the acceptance of the paper

Various other scripts are placed in "pipelines" folder.

### bargroupingIslands:
bargrouping_islands_synth (introduced version 0.3 - includes barcode islands plot)

calc_overlap_pcc_sort_m - does the number crunching of calculating overlaps between all pairs (OVERLAP)
create_barset - does the overlap graph layout (LAYOUT)
gen_assembly - plots barcode islands visually (if all barcodes have at least two clones) (CONSENSUS)
gen_reference_based_assembly - requires comparison to reference
plot_island - visually plots barcodes surrounding islands in different tabs (island have to have >=5 barcodes to be considered)

## real data:
bargrouping_islands_real_da (run for the chromosome assembly.)

### Explanation of parameters
1 | Plot figures | default.plotfigs
1 | Save figures | default.savefigs
0.95 | sFmin | default.sFmin
1.05 | sFmax | default.sFmax
0.01 | sFstep | default.sFstep
300 | minimum overlap | default.minOverlap
300 | minLen |  default.minLen
kymo | method |     default.method
20 | numframes |    default.numframes
0 | depth |    default.depth 
2 | Stouffer Z-score  cut-off| default.zscoreCutOff
hierarchical | comparison method | default.comparisonMethod
0.02 | Allowed scaling factor diff for merging | default.scDiffSetting
50 |    Allowed pixel difference for merging | default.pxDifSetting
[1 2]  | checklist items | clIdx
1 | whether use gui | useGUI
