function [F] = beta_ev_fast_mle_nu(x, cc)
    % Args:
    %   x - nu_eff
    %   cc - pcc values

    % Returns:
    %   F - score

    % Example


  % size of cc
  m = length(cc);

  % psi function difference
  psiDif = (psi((x - 1) / 2) - psi((x - 2) / 2));

  % now compute the second term

  term2 = m * psiDif;

  % finally the last term

  term3 = sum(log(1 - cc.^2));

  F = term2 + term3;
end
