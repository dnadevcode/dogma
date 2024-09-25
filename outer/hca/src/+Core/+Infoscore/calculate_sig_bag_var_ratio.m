function [bgmean,bgstd,sigbgstd,sigstd] = calculate_sig_bag_var_ratio(kymoStructs)

bgmean = zeros(1,length(kymoStructs));
bgstd = zeros(1,length(kymoStructs));
sigbgstd = zeros(1,length(kymoStructs));
sigstd = zeros(1,length(kymoStructs));

for i=1:length(kymoStructs)
    %background pixels
    kymo = double(kymoStructs{i}.alignedKymo);
    kymo(logical(kymoStructs{i}.alignedMask)) = nan;
    kymo(kymo==0) = nan;
    kymo = kymo(:);  
    kymo(isoutlier(kymo,'median')) = nan;
    kymo(kymo<0) = nan;

    % bg mean./std
    bgstd(i) = std(double(kymo(:)),1,'omitnan');
    bgmean(i) = mean(double(kymo(:)),'omitnan');

    % now signal
    kymo = double(kymoStructs{i}.alignedKymo);
    kymo(~logical(kymoStructs{i}.alignedMask)) = nan;
    kymo(kymo==0) = nan;
    kymo = kymo(:);
    kymo(isoutlier(kymo,'median')) = nan;

%     kymo(kymo < bgstd(i)+bgmean(i)) = nan; % remove outliers instead?

    sigbgstd(i) = std(double(kymo(:)),1,'omitnan');

    sigstd(i) = sqrt( sigbgstd(i).^2 -  bgstd(i) .^2);


end

