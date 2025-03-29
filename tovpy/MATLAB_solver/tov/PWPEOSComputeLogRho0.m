function lgrho0 = PWPEOSComputeLogRho0(G0,G1,K0,lgP1)

%PWPEOSComputeLogRho0 Compute log10(rho0[CGS]) from log10(P1[CGS])
%
%   lgrho0 = PWPEOSComputeLogRho0(G0,G1,K0,lgP1)
%
%   Assume log10(rho1[cgs]) = 14.7 
%   K0 must be given in c=G=Msun=1 units

Pcgs = 1.80171810e-39; 
Rcgs = 1.61930347e-18;
         
rho1 = 8.114744333450496e-04;  % -> lgrho[cgs] = 14.7                
K1   = (10^lgP1*Pcgs)/rho1^G1;                 
         
f    = @(x,K0,G1,K1)( K0*x.^G0 - K1*x.^G1 );      
lgrho0 = log10( fzero(@(x) f(x,K0,G1,K1), rho1 ) / Rcgs );