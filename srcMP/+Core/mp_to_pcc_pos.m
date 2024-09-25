function [cS] = mp_to_pcc_pos(cS)

for i=1:length(cS)
    cS{i}.pos =  cS{i}.pos- cS{i}.secondPos+1; % convert
end

