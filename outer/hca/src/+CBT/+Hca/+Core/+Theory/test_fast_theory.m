% later move this to TEST folder
import CBT.Hca.Core.Theory.compute_hca_theory_fast;

    sets = ini2struct( 'sets.txt');


chr = memmapfile('/home/albyback/git/hca/chr1.mm', 'format', 'uint8');   
 [ theoryCurveUnscaled_pxRes] = compute_hca_theory_fast( chr,sets);
