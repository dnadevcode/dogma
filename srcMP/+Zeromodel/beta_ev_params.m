function [parameters] = beta_ev_params(cc, x0, N)
  %compute_evd_params
  % Compute functional fomr parameters for max cc histogram
  %     Args:
  %       cc,x0
  %     Returns:
  %         parameters
  
  if nargin >=3
      f = @(y) abs(Zeromodel.beta_ev_fast_mle(y, cc, N));
  else
	f = @(y) abs(Zeromodel.beta_ev_fast_mle(y, cc));
  end

  % first parameter, the value should be found by running minimization
  % in the interval [2,x0]
  [x2,fval,exitflag,output] = fminbnd(f, 2, x0,optimset('TolX',1e-10));
% figure,plot(95:0.01:105,arrayfun(@(x) f(x),95:0.01:105))

  % second parameter - if N  provided, should be similar..
  m = length(cc);

    if nargin <3

      denom = (-1 / m * sum(log(1 + betainc(cc.^2, 1/2, x2 / 2 - 1))) + log(2));
      N2 = max(0, 1 ./ denom);
    else
        N2=N;
    end

  parameters = [x2 N2];
end
