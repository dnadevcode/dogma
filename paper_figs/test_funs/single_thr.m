ix=11
iy=3

 c = parcluster;
% c.AdditionalProperties.AccountName = 'snic2021-5-132';
% c.AdditionalProperties.WallTime = '1:00:00';

windowWidths = 350;

local_pipeline_mp(ix, iy, dirName, depth, windowWidths, sF, thryFiles)


[rezMaxMP] = hca_compare_distance(barGen(passingThreshBars), theoryStruct(12977), sets );
[rezMaxMP] = hca_compare_distance(barGen(passingThreshBars), theoryStruct(15846), sets );



