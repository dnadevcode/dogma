function [F] = beta_ev_fast_mle(nu, Cr, N)
  % beta_ev_fast_mle
  %     MLE function for beta ev using ContigAssembly formula S.M.43-S.M.44
  %     Defined for both Cr > 0 and Cr < 0

  %     Args:
  %         nu - nu_eff value
  %         cc - PCC coefficients
  %         N - N_eff

  %     Returns:
  %         F - equation for nu_eff = 0

  %     Example:
  %         import Zeromodel.beta_ev_fast_mle
  %         [F] = beta_ev_fast_mle(5, Cr, 100)

  % size of cc
  R = length(Cr);

  if nargin < 3
      % first calculate the value of N
      denom = (-1 / R * sum(log(1 + sign(Cr) .* betainc(Cr.^2, 1/2, nu / 2 - 1))) + log(2));
      N = 1 ./ denom;
  end

  %
  nVal = nu / 2 - 1;
  %nVal = x;

  % approximate the derivative here
  h = 0.00001;
  dVals = (log(1 + sign(Cr) .* betainc(Cr.^2, 1/2, nVal + h)) - log(1 + sign(Cr) .* betainc(Cr.^2, 1/2, nVal - h))) ./ (2 * h);

  term1 = (N - 1) * sum(dVals);

  % psi function difference
  psiDif = (psi((nu - 1) / 2) - psi((nu - 2) / 2));

  % now compute the second term

  term2 = R * psiDif;

  % finally the last term

  term3 = sum(log(1 - Cr.^2));

  F = term1 + term2 + term3;
end
