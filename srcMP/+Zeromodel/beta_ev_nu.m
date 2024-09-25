function [x2] = beta_ev_nu(cc, x0)
  %compute_evd_params
  % Compute functional fomr parameters for max cc histogram
  %     Args:
  %       cc,x0
  %     Returns:
  %         parameters

  f = @(y) abs(Zeromodel.beta_ev_fast_mle_nu(y, cc));


  % first parameter, the value should be found by running minimization
  % in the interval [2,x0]
  [x2,fval,exitflag,output] = fminbnd(f, 2, x0,optimset('TolX',1e-10));

%   % second parameter - if N  provided, should be similar..
%   m = length(cc);
% 
%     if nargin <3
% 
%       denom = (-1 / m * sum(log(1 + betainc(cc.^2, 1/2, x2 / 2 - 1))) + log(2));
%       N2 = max(0, 1 ./ denom);
%     else
%         N2=N;
%     end

%   parameters = [x2 N2];
end
