function [p,e,rho] = EOSGlaw(h, K,G)

%EOSGlaw gamma-law EOS (polytropic)
%
%  [p,e,rho] = eosGlaw(h, K,G)
%


Gmo = G-1;

H1 = 2.*sinh(0.5*h).*exp(0.5*h); % = exp(H) - 1;

rho = ( Gmo/G/K .*H1 ).^( 1.0/Gmo );

rhog = rho.^G;
p  = K*rho.^G;
e  = rho + p/Gmo;
