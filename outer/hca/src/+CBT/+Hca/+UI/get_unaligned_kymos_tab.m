function [hTabUnalignedKymos] = get_unaligned_kymos_tab(tsHCA)
    persistent localhTabUnalignedKymos;
    if isempty(localhTabUnalignedKymos) || not(isvalid(localhTabUnalignedKymos))
        hTabUnalignedKymos = tsHCA.create_tab('Unaligned Kymos');
        localhTabUnalignedKymos = hTabUnalignedKymos;
    else
        hTabUnalignedKymos = localhTabUnalignedKymos;
    end
end

