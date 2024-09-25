# DOGMA

requires SCAMP CLI for MP based bargrouping. (https://github.com/zpzim/SCAMP)

'''
git clone https://github.com/zpzim/SCAMP
cd SCAMP
git submodule update --init --recursive
mkdir build && cd build
'''

most recent version of HCA to compare alignment (https://gitlab.com/dnadev/hca)
most recent version of lldev for loading data into correct format (https://gitlab.com/dnadev/lldev)

The scripts are placed in "pipelines" folder.

### Running via GUI:
- bargi_gui
The default settings are in files/bargi_settings.txt

### Running via Script
- bargi_run / bargrouping_islands_real_da

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
