% Example
% f = @(y) Functions.rootFunction2( y, data );
            % x0 = [10];
             
 %            parameters = [fsolve(f,x0,optimoptions('fsolve','Display','off'))];
%



function [f] = full_dist_mle(cc, rho, n)
    % full distributions for CC values
    if rho ~= 0
%         f = (n-2)*gamma(n-1)*(1-rho.^2)^((n-1)/2)*(1-r.^2)^((n-4)/2))/...
%             (sqrt(2*pi)*gamma(n-1/2)*(1-rho*r)^(n-3/2))*hypergeom([1/2 1/2],(1/2)*(2*n-1),1/2*(rho*r+1));%
    else
        %f =  (1-cc^2)^((n-4)/2)/beta(1/2,1/2*(n-2));
        m = size(cc,1);
        C1 = 1/m*sum(log(1-cc.^2));        
        f = -psi((n-2)/2)+C1+psi((n-1)/2); 
    end



       
end

