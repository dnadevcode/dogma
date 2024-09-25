function [ F ] = functional_form_helping_fun(x, cc )
    % function for maximum likelihood stuff
    
    % size of the sample set
    % x is the value we test for, cc is the set of correlation coefficients
    
    % size of cc
    m = length(cc);
    
    % first calculate the value of N
    denom = (-1/m*sum(log(1+betainc(cc.^2,1/2,x/2-1)))+log(2));
    N = 1./denom;
    
    % 
    nVal = x/2-1;
    %nVal = x;
    
    % approximate the derivative here
	h = 0.00001;
	dVals= (log(1+betainc(cc.^2,1/2,nVal+h))-log(1+betainc(cc.^2,1/2,nVal-h)))./(2*h);
    
    term1 = (N-1)*sum(dVals);
    
    % psi function difference
    psiDif =( psi((x-1)/2)-psi((x-2)/2) );

    % now compute the second term
    
    term2 = m*psiDif;
    
    % finally the last term
    
    term3 =sum(log(1-cc.^2));
    
    F = term1+term2+term3; 
end

