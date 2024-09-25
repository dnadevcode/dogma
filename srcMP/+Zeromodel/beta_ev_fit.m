function [a_fit, n_fit] = beta_ev_fit(y, lowerBound, upperBound, startParams, isConstant)

  if nargin < 2
    lowerBound = [4 1];
  end

  if nargin < 3
    upperBound = [inf inf];
  end

  if nargin < 4
    startParams = min(max([4, 1], lowerBound), upperBound);
  end

  if nargin < 5
    isConstant = false(1, 2);
  end

  import Zeromodel.beta_ev_pdf
  import Zeromodel.beta_ev_cdf

  switch bin2dec(num2str(isConstant))
    case 0
      pdfit = @(data, a, n) beta_ev_pdf(data, a, 1, n);
      cdfit = @(data, a, n) beta_ev_cdf(data, a, 1, n);
    case 1
      pdfit = @(data, a) beta_ev_pdf(data, a, 1, startParams(2));
      cdfit = @(data, a) beta_ev_cdf(data, a, 1, startParams(2));
    case 2
      pdfit = @(data, n) beta_ev_pdf(data, startParams(1), 1, n);
      cdfit = @(data, n) beta_ev_cdf(data, startParams(1), 1, n);
    case 3
      error("At least one parameter must not be constant.")
  end

  warning('off', 'stats:mlecov:NonPosDefHessian')
  warning('off', 'stats:mle:IterLimit')
  opt = statset('MaxIter', 1e6, 'MaxFunEvals', 1e6, 'TolX', 1e-6); %,'Display', 'iter');
  [params, ~] = mle(y(~isnan(y)), 'pdf', pdfit, ...
    'start', startParams(not(isConstant)), ...
    'cdf', cdfit, ...
    'Lowerbound', lowerBound(not(isConstant)), ...
    'Upperbound', upperBound(not(isConstant)), ...
    'Alpha', 0.01, ...
    'Options', opt);
  warning('on', 'stats:mlecov:NonPosDefHessian')
  warning('on', 'stats:mle:IterLimit')

  a_fit = startParams(1);
  n_fit = startParams(2);

  switch bin2dec(num2str(isConstant))
    case 0
      a_fit = params(1);
      n_fit = params(2);
    case 1
      a_fit = params(1);
    case 2
      n_fit = params(1);
  end
