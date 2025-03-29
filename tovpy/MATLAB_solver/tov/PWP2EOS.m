function eos = PWP2EOS(name)

%PWP2EOS Parameters of 2-pieces piecewise polytrope EOS
%
%   eos = PWP2EOS(name)
%
%   Refs
%   T.Damour & A.Nagar, Phys.Rev.D80:084035,2009
%   Uryu et al, Phys.Rev.D80:124004,2009
%   J.S.Read et al Phys.Rev.D79:124032,2009
%   
%   Adapted from PWP4EOS.m to account for results in
%   Lackey et al. 1109.3402
%   see Tab. I there
%

%  Notes
%     
% 1. all the PWP EOS uses
% 
%    K0 = 3.5966*1.e13;  % [cgs]
%    G0 = 1.3569;
%
%    if others are needed, modify inside the case distinction
%
% 2. diving log densities 14.7 and 15 are typically used, corresponding to
%
%      rhonuc = 2.7e14 , log10( [1.85 3.70]*rhonuc )
%
% 3. rhoc refers to a particular test model with mass = 1.3500(0)
%
% 4. rhom refers to the maximum mass model
%
% 5. Usually the log10(P1[cgs]) is given instead of lgrho0 (used here)
%
%     you compute lgrho0 [cgs] like this, e.g.
%
%         lgP1 = 34.269;
%
%         Pcgs = 1.80171810e-39; Rcgs = 1.61930347e-18;
%         rho1 = 8.114744333450496e-04;  % -> lgrho[cgs] = 14.7
%         K1   = (10^lgP1*Pcgs)/rho1^G1;
%         K0   = 8.949298583670957e-02;  % -> K0 [cgs] = 3.59389*1.e13
%         f    = @(x,K0,G1,K1)( K0*x.^G0 - K1*x.^G1 );
%         log10( fzero(@(x) f(x,K0,G1,K1), rho1 ) / Rcgs )
%
%    in this version a routine takes care of this.


% Unit conversion
GNewt     = 6.67384*1e-8;     % cm^3 g^-1 s^-2
c         = 2.99792458*1e+10; % cm/sec
c2        = c*c;
GMoc2     = 0.5*2.953250077*1e5; % 1.476701332464468e+05
GMoc3     = GMoc2/c;
Msun      = GMoc2;

lgc = 0.00364; % correction term for lgP1, see Lackey et al. 1109.3402 Tab. I


% Constants
G0 = 1.3569; % <--------------------------------------------- G0
cgs.K0 = 3.5966*1e13; % [cgs] % <------------------------- K0_cgs


% Determine K in units of length: [K] = [L]^{2(G0-1)}
K  = GNewt/(c^4)*( c^2/GNewt )^G0*cgs.K0;
K0 = K/(Msun^(2*(G0-1)));

rho0 = []; lgrho0 = [];


% EOS specific pars
switch name
    
    case 'p.3G2.4_Bss'
        
        lgP1 = 34.3;
        G1 = 2.4;
        
        rho0 = 1.383616e-04;
        lgrho0 = 1.393168730414e+01;
        
        lgrhoc = 1.517450940656e+01;
        lgrhom = 1.546745400738e+01;
        
    case 'p.3G2.7_Bs'
        
        lgP1 = 34.3;
        G1 = 2.7;
        
        rho0 = 2.054078e-04;
        lgrho0 = 1.410328866457e+01;
        
        lgrhoc = 1.505643802511e+01;
        lgrhom = 1.539659815011e+01;
        
    case 'p.3G3.0_B'
        
        lgP1 = 34.3;
        G1 = 3.0;
        
        rho0 = 2.639707e-04;
        lgrho0 = 1.421222748696e+01;
        
        lgrhoc = 1.499013997254e+01;
        lgrhom = 1.533782014956e+01;
        
    case 'p.3G3.3'
        
        lgP1 = 34.3;
        G1 = 3.3;
        
        rho0 = 3.139467e-04;
        lgrho0 = 1.428752764256e+01;
        
        lgrhoc = 1.494613665483e+01;
        lgrhom = 1.528430089383e+01;
        
    case 'p.4G2.4_HBss'
        
        lgP1 = 34.4;
        G1 = 2.4;
        
        rho0 = 1.109552e-04;
        lgrho0 = 1.383581921862e+01;
        
        lgrhoc = 1.503503692979e+01;
        lgrhom = 1.539106566151e+01;
        
    case 'p.4G2.7_HBs'
        
        lgP1 = 34.4;
        G1 = 2.7;
        
        rho0 = 1.730463e-04;
        lgrho0 = 1.402883404466e+01;
        
        lgrhoc = 1.496220071264e+01;
        lgrhom = 1.533782014956e+01;
        
    case 'p.4G3.0_HB'
        
        lgP1 = 34.4;
        G1 = 3.0;
        
        rho0 = 2.294538e-04;
        lgrho0 = 1.415136691852e+01;
        
        lgrhoc = 1.491777264478e+01;
        lgrhom = 1.528784617460e+01;
        
    case 'p.4G3.3'
        
        lgP1 = 34.4;
        G1 = 3.3;
        
        rho0 = 2.788635e-04;
        lgrho0 = 1.423606348734e+01;
        
        lgrhoc = 1.488710790462e+01;
        lgrhom = 1.524326404195e+01;
        
    case 'p.5G2.4'
        
        lgP1 = 34.5;
        G1 = 2.4;
        
        rho0 = 8.897733e-05;
        lgrho0 = 1.373995113311e+01;
        
        lgrhoc = 1.490738915538e+01;
        lgrhom = 1.532504586152e+01;
        
    case 'p.5G2.7'
        
        lgP1 = 34.5;
        G1 = 2.7;
        
        rho0 = 1.457833e-04;
        lgrho0 = 1.395437942475e+01;
        
        lgrhoc = 1.487134470387e+01;
        lgrhom = 1.527712230921e+01;
        
    case 'p.5G3.0_H'
        
        lgP1 = 34.5;
        G1 = 3.0;
        
        rho0 = 1.994503e-04;
        lgrho0 = 1.409050635008e+01;
        
        lgrhoc = 1.484685223674e+01;
        lgrhom = 1.523933374832e+01;
        
    case 'p.5G3.3'
        
        lgP1 = 34.5;
        G1 = 3.3;
        
        rho0 = 2.477009e-04;
        lgrho0 = 1.418459933213e+01;
        
        lgrhoc = 1.482883353119e+01;
        lgrhom = 1.519794106316e+01;
        
%     case 'p.6G2.4'
%         
%         lgP1 = 34.6;
%         G1 = 2.4;
%         
%         %rho0 = 8.936925e-04;
%         %lgrho0 = 1.474185987755e+01;
%         
%         lgrhoc = 15.00; % fixme
%         lgrhom = 15.00; 
%         
    case 'p.6G2.7'
        
        lgP1 = 34.6;
        G1 = 2.7;
        
        rho0 = 1.228155e-04;
        lgrho0 = 1.387992480484e+01;
        
        lgrhoc = 1.478301151668e+01; 
        lgrhom = 1.522324692842e+01;
        
    case 'p.6G3.0'
        
        lgP1 = 34.6;
        G1 = 3.0;
        
        rho0 = 1.733700e-04;
        lgrho0 = 1.402964578164e+01;
        
        lgrhoc = 1.477711172539e+01;
        lgrhom = 1.518471279743e+01;
        
    case 'p.6G3.3'
        
        lgP1 = 34.6;
        G1 = 3.3;
        
        rho0 = 2.200206e-04;
        lgrho0 = 1.413313517691e+01;
        
        lgrhoc = 1.477120439204e+01; 
        lgrhom = 1.515218357260e+01;
        
%     case 'p.7G2.4'
%         
%         lgP1 = 34.7;
%         G1 = 2.4;
%         
%         %rho0 = 7.166713e-04;
%         %lgrho0 = 1.464599179204e+01;
%                 
%         lgrhoc = 15.00; % fixme
%         lgrhom = 15.00;
        
    case 'p.7G2.7'
        
        lgP1 = 34.7;
        G1 = 2.7;
        
        rho0 = 1.034662e-04;
        lgrho0 = 1.380547018493e+01;
        
        lgrhoc = 1.469664081487e+01;
        lgrhom = 1.516172889050e+01;
        
    case 'p.7G3.0_1.5H'
        
        lgP1 = 34.7;
        G1 = 3.0;
        
        rho0 = 1.507001e-04;
        lgrho0 = 1.396878521321e+01;
        
        lgrhoc = 1.470835383750e+01;
        lgrhom = 1.513746031578e+01;
        
    case 'p.7G3.3'
        
        lgP1 = 34.7;
        G1 = 3.3;
        
        rho0 = 1.954336e-04;
        lgrho0 = 1.408167102169e+01;
        
        lgrhoc = 1.471413429543e+01;
        lgrhom = 1.511175491554e+01;
        
    case 'p.9G3.0_2H'
        
        lgP1 = 34.9;
        G1 = 3.0;
        
        rho0 = 1.138655e-04;
        lgrho0 = 1.384706407633e+01;
        
        lgrhoc = 1.457321205829e+01;
        lgrhom = 1.503679015390e+01;
        
    otherwise
        
        % add more from Tab III of
        % http://arxiv.org/pdf/0812.2163.pdf
        
        error('unknown EOS %s',name)
end


% Dimensionless rho's
lgrho2  = 15.0;
lgrho1  = 14.7;

cgs.rho2 = 10^(lgrho2); % g/cm3
cgs.rho1 = 10^(lgrho1); % g/cm3

rho2 = cgs.rho2 *GNewt/c2*Msun.^2;
rho1 = cgs.rho1 *GNewt/c2*Msun.^2;

if isempty(rho0)
    % compute the rho0 is not pre-computed (rho1 is a guess)
    [rho0,cgs.rho0,lgrho0] = lgP1_to_lgrho0( lgP1 + lgc, K0,rho1, G0,G1 );
    %[rho0,cgs.rho0,lgrho0] = lgP1_to_lgrho0( lgP1 + lgc, K0,rho1*1.3, G0,G1 );
    fprintf(' * %s\n rho0 = %.6e;\n lgrho0 = %.12e;\n',name,rho0, lgrho0);
end

cgs.rho0 = 10^(lgrho0);

cgs.rhoc = 10^(lgrhoc); % g/cm3
cgs.rhom = 10^(lgrhom); % g/cm3

rhoc = cgs.rhoc *GNewt/c2*Msun.^2;
rhom = cgs.rhom *GNewt/c2*Msun.^2;


% Other quantities
K1 = K0*rho0^(G0-G1);

a1 = K0*rho0^(G0-1)*( 1/(G0-1) - 1/(G1-1) );
e1 = (1+a1)*rho1 + K1/(G1-1)*rho1.^G1;
a2 = e1/rho1 - 1 - K1/(G1-1)*rho1.^(G1-1);
e2 = (1+a2)*rho2 + K1/(G1-1)*rho2.^G1;
a3 = e2/rho2 - 1 - K1/(G1-1)*rho2.^(G1-1);


% Output
% put in pwp4 form as this is what the TOVL.m can read
eos.type = 'pwp4';
eos.name = name;
eos.rho  = [rho0 rho1 rho2];
eos.G    = [G0 G1 G1 G1];
eos.K    = [K0 K1 K1 K1];
eos.a    = [a1 a2 a3];
eos.e    = [e1 e2   ];
eos.rhoc = rhoc;
eos.rhom = rhom;




function [rho,rho_cgs,lgrho0] = lgP1_to_lgrho0( lgP1, K0,rho1, G0,G1 )

% Routine to convert from lgP1 to lgrho0

Pcgs = 1.80171810e-39;
Rcgs = 1.61930347e-18;
%rho1 = 8.114744333450496e-04; % -> lgrho[cgs] = 14.7 

K1   = (10^lgP1*Pcgs)/rho1^G1;

f    = @(x,K0,G1,K1)( K0*x.^G0 - K1*x.^G1 );
rho  = fzero(@(x) f(x,K0,G1,K1), rho1 );

rho_cgs = rho / Rcgs;
lgrho0 = log10( rho_cgs );

