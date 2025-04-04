function varargout = TOVh(hc, npts, eos)

%TOVh Simple TOV solver using Lindbloom formulation with enthalpy
% 
%   [m,r,x,phi,rho,p,h,mstar,rstar,mbar,R] = TOVh(hc, npts, eos)
%
%      hc    central enthalpy
%      npts  no radial points
%      eos   function handle to EOS routine
%
%            [p,e,rho] = eos(h)
%
%      m        mass function
%      r,x      Schwarzschild and isotropic radii
%      phi      other metric function
%      rho,p,h  fluid quantities 
%      mstar    star mass
%      rstar    star radius
%      mbar     rest-mass 
%      R        proper radius
%
%   sol = tov(hc, npts, eos) return a structure
%
%   Example
%          sol = TOVh(0.23, 1000, @(x) eosGlaw(x,100,2)); 
%

% Spacing
dh = hc/npts;

% Central values from EoS
[pc,ec] = eos(hc);

% Initial data
h     = hc - dh;        
[p,e] = eos(h);
decdh = (e-ec)/(h-hc);

r = sqrt(3.*dh/(2.*pi*(ec+3.*pc)))* ...
    (1.-0.25*(ec-3.*pc-3./5.*decdh)*(dh/(ec+3.*pc)));

m = 4./3.*pi*ec*r*r*r*(1.-3./5.*decdh*dh/ec);

y0 = [r, m, 0];

% h-grid
hgrid = [h:-dh:0];

% Matter fields
[p,e,rho] = eos(hgrid);

% Integrate
options = odeset('RelTol',1e-9,'AbsTol',1e-12);
[h,y] = ode45(@(x,y) toveqs(x,y, eos),hgrid,y0,options);

% radius and mass
r  = y(:,1);
m  = y(:,2);
Ir = y(:,3);

IR    = y(end,3);
rstar = r(end);
mstar = m(end);
phiR  = 0.5*log(1.-2.*mstar/rstar);

% Metric function
phi = phiR - h;

% Baryonic mass
mbar = barmass(r,m,rho');

% Proper radius
R = proprad(r,m);

% isotropic radius
C = 1./(2*r) .* (sqrt(r.*r-2*m.*r)+r-m) * exp(-IR);
x = r .* C .* Ir;

fprintf(1,' mass = %e\n radius = %e\n Mb = %e\n R = %e\n\n',...
    mstar,rstar,mbar,R);

if nargout>1
    varargout = {m,r,x,phi,rho,p,hgrid,mstar,rstar,mbar,R};
else
    TOV.m   = m;
    TOV.r   = r;
    TOV.x   = x;
    TOV.phi = phi;
    TOV.rho = rho;
    TOV.p   = p;
    TOV.h   = h;
    TOV.mstar = mstar;
    TOV.rstar = rstar;
    TOV.mbar = mbar;
    TOV.R   = R;
    varargout = {TOV};
end





function dy = toveqs(h,y, eos)

dy = zeros(3,1);

[p,e] = eos(h);

r    = y(1);
m    = y(2);
r3   = r*r*r;
f    = sqrt(1-2*m/r);

tmp  = - (r-2.*m)/(m+4*pi*r3*p);
drdh = r *tmp;
dmdh = 4*pi*e*r3 *tmp;
dIdh = (1-f)/(r*f)*drdh;

dy = [drdh; dmdh; dIdh];
  


function mb = barmass(r,m,rho)
tmp = r.^2.*rho./sqrt(1-2*m./r);
mb = 4*pi*trapz(r,tmp);



function pr = proprad(r,m)
pr = trapz(r, 1./sqrt(1-2*m./r));




%{ 
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
%}


