function [kymoStructs] = kymo_sig_std(kymoStructs)
    % kymo_sig_std - estimates kymo signal and background standard
    % deviation

    for idFold=1:length(kymoStructs)
        for idKym=1:length(kymoStructs{idFold})
            kymo = double(kymoStructs{idFold}{idKym}.unalignedKymo);
            kymo(logical(kymoStructs{idFold}{idKym}.unalignedBitmask)) = nan;
            kymo(kymo==0) = nan;
            kymo(isoutlier(kymo)) = nan;
            % succeptible to outlier calculation
            kymoStructs{idFold}{idKym}.bgstd = std(double(kymo(:)),1,'omitnan');
            kymo = double(kymoStructs{idFold}{idKym}.unalignedKymo);
            kymo(logical(~kymoStructs{idFold}{idKym}.unalignedBitmask)) = nan;
            kymo(kymo==0) = nan;
            kymo(isoutlier(kymo)) = nan;
            kymo(kymo > mean(kymo,'omitnan')+std(double(kymo(:)),1,'omitnan')) = nan;
            kymoStructs{idFold}{idKym}.sigbgstd = std(double(kymo(:)),1,'omitnan');

            kymoStructs{idFold}{idKym}.sigstd = sqrt( kymoStructs{idFold}{idKym}.sigbgstd.^2- kymoStructs{idFold}{idKym}.bgstd.^2);
            if kymoStructs{idFold}{idKym}.sigbgstd < kymoStructs{idFold}{idKym}.bgstd
                kymoStructs{idFold}{idKym}.sigstd = nan;
            end
        end
    end

end

