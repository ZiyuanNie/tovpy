function k = ComputeLegendre(c,y,l)
% ComputeLegendre computes the Love number kl
%   Given the compactness c, y = R H(R)'/H(R) and the multipole index
%   l, compute Love number kl.

% Convert c to x.
x = 1/c - 1;

% Compute normalization coefficients for Legendre functions.
nP = -prod((2*l-1)/2-(0:(l-1)))/factorial(l) * 2^l * l*(l-1);
nQ = factorial(l)/prod(2*(0:l)+1)*(l+1)*(l+2);

% Compute Legendre function of the first kind and its derivative with respect
% to x.
Pl2  = 0;
Pl2p = 0;
for i=2:l
    Pl2  = Pl2  + factorial(i)/factorial(i-2) * nchoosek(l,i) * prod((l+i-1)/2-(0:(l-1)))/factorial(l) * x^(i-2);
    Pl2p = Pl2p + factorial(i)/factorial(i-2) * nchoosek(l,i) * prod((l+i-1)/2-(0:(l-1)))/factorial(l) * (i-2)*x^(i-3); 
end
Pl2p = 2^l*(-2*x)*Pl2/nP + 2^l*(1-x^2)*Pl2p/nP;
Pl2  = 2^l*(1-x^2)*Pl2/nP;

% Compute Legendre function of the second kind and its derivative with
% respect to x.
Ql2  = 1/nQ * sqrt(pi)/2^(l+1) * gamma(l+3)/gamma(l+3/2) * (x^2-1)/x^(l+3) * hypergeom([(l+3)/2 (l+4)/2],l+3/2,1/x^2);
Ql2p = 1/nQ * sqrt(pi)/2^(l+1) * gamma(l+3)/gamma(l+3/2) * (2*x^(-2 - l)*hypergeom([(l+3)/2 (l+4)/2],l+3/2,1/x^2) +...
                                                            (-3 - l)*x^(-4 - l)*(-1 + x^2)*hypergeom([(l+3)/2 (l+4)/2],l+3/2,1/x^2) -...
                                                            (2*((l+3)/2)*((l+4)/2)*x^(-6 - l)*(-1 + x^2)*hypergeom([(l+3)/2+1 (l+4)/2+1],l+3/2+1,1/x^2)/(l+3/2)));
                                                        
% Compute Love number kl.
k = -1/2 * c^(2*l+1) * (Pl2p - c*y*Pl2) / (Ql2p - c*y*Ql2);

end

