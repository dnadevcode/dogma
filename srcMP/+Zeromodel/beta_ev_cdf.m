function [cdfVal] = beta_ev_cdf(x, a, b, n, extraPrecision)
    % beta_ev_cdf
    % Extreme value cummulative distribution function for generalized beta

    %   Args:
    %   x - CC value
    %   a - parameter a   (nu_eff )
    %   b - parameter b    (typically 1)
    %   n - number of placement attempts

    %   Returns:
    %       cdfVal - cdfVal value at x

    %   Example:
    %     import Zeromodel.beta_ev_cdf;
    %     [cdfVal] = arrayfun(@(x) beta_ev_cdf(x, 80, 1, 1000),0.4:0.01:0.9);
    %     figure,plot(0:0.01:1,cdfVal);
    %     xlabel('CC value')
    %     ylavel('CDF')

    undefx = x.^2 <= 0 | x.^2 >= 1 | isnan(x);
    x(undefx) = 0.5; % note: should skip undefined values alltogether
    
    if nargin < 5 ||~extraPrecision
        cdfVal = (0.5 * (1 + sign(x) .* betainc(x.^2, 0.5 * b, a / 2 - 1, 'lower'))).^n;
    else
        cdfVal = (0.5 * (vpa(2) - sign(x) .* betainc(x.^2, 0.5 * b, a / 2 - 1, 'upper'))).^n;
        cdfVal = min(cdfVal, vpa(1) - 10^ - digits);
    end
    
    cdfVal(undefx) = nan;
    
    if nargin > 4 && extraPrecision
        cdfVal(x == 0) = 10^ - digits;
        cdfVal(x == 1) = vpa(1) - 10^ - digits;
    else
        cdfVal(x == 0) = 0;
        cdfVal(x == 1) = 1;
    end
end