function [] = save_disc(barcodeGenC,thryNames,refNums,rezMax,bestCoefs,is_distinct,refNumBad,bestCoefsBad,hcaSets)


plotNonDistinct = 1;
outFold = fullfile(fileparts(fileparts(hcaSets.kymofolder{1})),['distinct',hcaSets.timestamp,'.csv']);

  fid = fopen(outFold, 'w');
  fprintf(fid, '"Barcode_name","First_match", "First_match_CC","len_first_match", "Second_match", "Second_match_CC","len_second_match","Is_distinct","Match_ref_idx","Match_ref2_idx"\n');
  for i=1:length(is_distinct)
      if is_distinct(i)==1 || plotNonDistinct
          [~,fname] = fileparts(barcodeGenC{i}.name);
          first_match = strrep(thryNames{refNums{i}(1)},',','');
          first_match_cc = bestCoefs{i}(1);
          len_first_match = rezMax{refNums{i}(1)}{i}.lengthMatch;
          if~isnan(refNumBad)
              second_match = strrep(thryNames{refNumBad(i)},',','');
              second_match_cc = bestCoefsBad(i);
                len_second_match = rezMax{refNumBad(i)}{i}.lengthMatch;
          else
            second_match = '';
            second_match_cc = nan;
            len_second_match = nan;
          end

          fprintf(fid,'%s,%s,%4.3f,%4.3f,%s,%4.3f,%4.3f,%4.3f,%s,%4.3f\n',fname,first_match,first_match_cc,len_first_match,second_match,second_match_cc,len_second_match,is_distinct(i),num2str(refNums{i}),refNumBad(i));
      end
  end
  fclose(fid);
    display(['Saved output of discriminatory analysis to ',outFold]);


end

