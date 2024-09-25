
function tableSynth = synth_struct_to_table(synthStr)
    tableSynth = zeros(length(synthStr),4);
    for i=1:length(synthStr)
        tableSynth(i,:) = [synthStr{i}.pos synthStr{i}.pos+synthStr{i}.lengthMatch-1 synthStr{i}.or synthStr{i}.rf];
    end

end