function [f] = full_dist(cc, rho, n)

   if rho ~= 0
%         f = (n-2)*gamma(n-1)*(1-rho.^2)^((n-1)/2)*(1-r.^2)^((n-4)/2))/...
%             (sqrt(2*pi)*gamma(n-1/2)*(1-rho*r)^(n-3/2))*hypergeom([1/2 1/2],(1/2)*(2*n-1),1/2*(rho*r+1));%
    else
       f =  (1-cc^2)^((n-4)/2)/beta(1/2,1/2*(n-2));
   end

end

