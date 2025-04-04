function [Lamt] = LambdaTilde(Lam1,Lam2, M1,M2)

%LambdaTilde Compute the beloved tidal parameter \tilde{\Lambda} 
%
%   [Lamt] = LambdaTilde(Lambda1,Lmabda2, M1,M2)
%

Lamt = 16/13*( ...
    (M1 + 12*M2)/(M1 + M2)^5 * M1^4 * Lam1 + ...
    (M2 + 12*M1)/(M1 + M2)^5 * M2^4 * Lam2   ...
    );
