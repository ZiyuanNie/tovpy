function [Lam] = Lambda(kl,C,l)

%Lambda Compute tidal parameter Lambda from Love numbers
%
%   Lam = Lambda(kl,C, l)
%   k(l) multipolar Love numbers, l=1...l_max 
%   C compactness 
%   (k(l) and l must be same size!)
%

f2lm1 = Factd(2*l-1);
Cpow2l1 = C.^(2*l+1);
Lam = 2 * kl ./ (Cpow2l1.*f2lm1);

Lam(l<2) = 0; % dummy
