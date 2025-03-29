function eos = PWP4EOS(name)

%PWP4EOS Parameters of 4-pieces piecewise polytrope EOS
%
%   eos = PWP4EOS(name)
%
%   Refs
%   T.Damour & A.Nagar, Phys.Rev.D80:084035,2009
%   Uryu et al, Phys.Rev.D80:124004,2009
%   J.S.Read et al Phys.Rev.D79:124032,2009
%

%  Notes
%     
% 1. all the PWP EOS uses
% 
%    K0 = 3.59389*1.e13;  % [cgs]
%    G0 = 1.35692;
% 
%    if others are needed, modify inside the case distinction
%     
% 2. diving log densities 14.7 and 15 are typically used, corresponding to 
%  
%      rhonuc=2.7e14 , log10( [1.85 3.70]*rhonuc )
% 
% 3. rhoc refers to a particular test model with mass = 1.3500(0)
%  
% 4. rhom refers to the maximum mass model 
%  
% 5. Usually the log10(P1[cgs]) is given instead of lgrho0 (used here)
% 
%    you compute lgrho0 [cgs] like this:
%
%         % define G0, G1, G2, K0 and lgP1
% 
%         K0   = 8.949298583670957e-02;  % -> K0 [cgs] = 3.59389*1.e13
%         lgP1 = 34.269; 
% 
%         Pcgs = 1.80171810e-39; Rcgs = 1.61930347e-18;
%         rho1 = 8.114744333450496e-04;  % -> lgrho[cgs] = 14.7       
%         K1   = (10^lgP1*Pcgs)/rho1^G1;        
%         f    = @(x,K0,G1,K1)( K0*x.^G0 - K1*x.^G1 );
%         log10( fzero(@(x) f(x,K0,G1,K1), rho1 ) / Rcgs )        
% 
%    see PWPEOSComputeLogRho0.m

G0 = 1.35692; % <--------------------------------------------- G0
cgs.K0 = 3.59389*1e13; % [cgs] % <------------------------- K0_cgs
% K0 = 8.949298583670957e-02

% all EOS implemented:
eosname = {'ALF2','BGN1H1','ENG','FPS','H4','MPA1','MS1', ...
	   'MS1B','2B','2H','HB','SLY','APR4','H3','PS'};

% EOS specific pars
switch upper(name)
    
    case 'ALF2'
        
        G1 = 4.070;
        G2 = 2.411;
        G3 = 1.890;
                
        lgrho0 = 14.28950292;
        lgrho1 = 14.7;
        lgrho2 = 15.0;

        
        lgrhoc = 1.4805353669347676e+01;
        lgrhom = 1.521189490426865e+01; % max = 1.9909e+00
        
    case 'BGN1H1'
        
        G1 = 3.258;
        G2 = 1.472;
        G3 = 2.464;
        
        lgrho0 = 14.1110;
        lgrho1 = 14.7;
        lgrho2 = 15.0;

        
        lgrhoc = 1.4915206642657735e+01;
        lgrhom = 1.540656939097902e+01; % max = 1.6437e+00

    case 'ENG'
        
        G1 = 3.514;
        G2 = 3.130;
        G3 = 3.168;
        
        lgrho0 = 14.26667744;
        lgrho1 = 14.7;
        lgrho2 = 15.0;

        
        lgrhoc = 1.4882151441114733e+01;
        lgrhom = 1.524857545282772e+01; % max = 2.2505e+00
                       
    case 'FPS'
        
        G1 = 2.985;
        G2 = 2.863;
        G3 = 2.600;
        
        lgrho0 = 14.220;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.5036795442475308e+01;
        lgrhom = 1.537921917778453e+01; % max = 1.8101e+00
        
    case 'H4'
        
        G1 = 2.909;
        G2 = 2.246;
        G3 = 2.144;
                
        lgrho0 = 13.94829169;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.4741778853460449e+01;
        lgrhom = 1.520417790584642e+01; % max = 2.0281e+00
                
    case 'MPA1'
        
        G1 = 3.446;
        G2 = 3.572;
        G3 = 2.887;
                
        lgrho0 = 14.22480928;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.4820125068314582e+01;
        lgrhom = 1.517230864091296e+01; % max = 2.4703e+00
        
    case 'MS1'
        
        G1 = 3.224;
        G2 = 3.033;
        G3 = 1.325;
                
        lgrho0 = 13.9738868;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.4622051722343318e+01;
        lgrhom = 1.505341654098932e+01; % max = 2.7693e+00
        
    case 'MS1B'
        
        G1 = 3.456;
        G2 = 3.011;
        G3 = 1.425;
        
        lgrho0 = 14.05556938;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        lgrhoc = 1.4637283615473175e+01; 
        lgrhom = 1.506053054389700e+01; % max = 2.7628e+00
        
    case '2B'
        
        G1 = 3.0;
        G2 = 3.0;
        G3 = 3.0;
                
        lgrho0 = 14.334;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.5140917542011509e+01; 
        lgrhom = 1.543648569646101e+01; % max = 1.7833e+00
        
    case '2H'
        
        G1 = 3.0;
        G2 = 3.0;
        G3 = 3.0;
                
        lgrho0 = 13.847;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.4573115091911582e+01; 
        lgrhom = 1.504293816237676e+01; % max = 2.8345e+00
        
    case 'HB'
        
        G1 = 3.0;
        G2 = 3.0;
        G3 = 3.0;
                
        lgrho0 = 14.;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.4743477013547377e+01; 
        lgrhom = 1.516260737577637e+01; % max = 2.4510e+00
        
    case 'SLY'
        
        G1 = 3.005;
        G2 = 2.988;
        G3 = 2.851;
                
        lgrho0 = 14.165;
        lgrho1 = 14.7;
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.4933483482249795e+01; 
        lgrhom = 1.530455239900855e+01; % max = 2.0604e+00
            
    case 'APR4'
        
        G1 = 2.830;
        G2 = 3.445;
        G3 = 3.348;                 
         
        lgrho0 = 14.1790;
        lgrho1 = 14.7; 
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.4948397261990449e+01; 
        lgrhom = 1.528083677639307e+01; % max = 2.2006e+00

    case 'H3'
        
        G1 = 2.787;
        G2 = 1.951;
        G3 = 1.901;
                
        lgrho0 = 13.9276;
        lgrho1 = 14.6985; 
        lgrho2 = 14.9996;
        
        
        lgrhoc = 1.4834594319754567e+01; 
        lgrhom = 1.525784124817839e+01; % max = 1.7010e+00

    case 'PS'
        
        G1 = 2.216;
        G2 = 1.640;
        G3 = 2.365;
                
        lgrho0 = 13.907;
        lgrho1 = 14.7; 
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.2e+01; % ?
        lgrhom = 1.525570337100e+01; % ?   

  case 'KIUCHI19_3765'

        % https://arxiv.org/pdf/1903.01466.pdf 
        % 3 segment EOS with varying G1 and lgP1
        % Needs to redefine the 0 values
        G0 = 1.357;
        K0 = 8.951e-2;
        cgs.K0 = 3.582812504842164e+13; % = K0 *Msun^(2*(G0-1)) /( G/(c^4)*( c^2/G )^G0 )
        
        G1 = 3.765; 
        G2 = 2.8;
        G3 = 2.8;
                
        lgrho0 = 1.442648953305288e+01; % = PWPEOSComputeLogRho0(1.375,3.765,8.951e-2,34.1)
        lgrho1 = 14.7; 
        lgrho2 = 15.0;
        
        
        lgrhoc = 1.515776986084911e+01; 
        lgrhom = 1.545782295587499e+01; % this is 1.7Mo , does not match the 2Mo of the paper! Radius 1.35 also does not match

    case 'ALL'

        eos = eosname;
        return;

    otherwise
        
        % add more from Tab III of
        % http://arxiv.org/pdf/0812.2163.pdf
        
        error('unknown EOS %s',name)
        
end


% Unit conversion
G         = 6.67384*1e-8;     % cm^3 g^-1 s^-2
c         = 2.99792458*1e+10; % cm/sec
c2        = c*c;
GMoc2     = 0.5*2.953250077*1e5; % 1.476701332464468e+05
GMoc3     = GMoc2/c;
Msun      = GMoc2;

cgs.rho0 = 10^(lgrho0); % g/cm3
cgs.rho1 = 10^(lgrho1); % g/cm3
cgs.rho2 = 10^(lgrho2); % g/cm3

cgs.rhoc = 10^(lgrhoc); % g/cm3
cgs.rhom = 10^(lgrhom); % g/cm3


% Determine K in units of length: [K] = [L]^{2(G0-1)}
K  = G/(c^4)*( c^2/G )^G0*cgs.K0;
K0 = K/(Msun^(2*(G0-1)));


% Dimensionless rho's
rho0 = cgs.rho0*G/c2*Msun.^2;
rho1 = cgs.rho1*G/c2*Msun.^2;
rho2 = cgs.rho2*G/c2*Msun.^2;

rhoc = cgs.rhoc*G/c2*Msun.^2; 
rhom = cgs.rhom*G/c2*Msun.^2; 


% Other quantities
K1 = K0*rho0^(G0-G1);
K2 = K1*rho1^(G1-G2);
K3 = K2*rho2^(G2-G3);

a1 = K0*rho0^(G0-1).*( 1/(G0-1) - 1/(G1-1) );
e1 = (1+a1)*rho1 + K1/(G1-1).*rho1.^G1;
a2 = e1/rho1 - 1 - K1/(G2-1)*rho1.^(G1-1);
e2 = (1+a2)*rho2 + K2/(G2-1).*rho2.^G2;
a3 = e2/rho2 - 1 - K2/(G3-1)*rho2.^(G2-1);


% Output
eos.type = 'pwp4';
eos.name = name;
eos.rho  = [rho0 rho1 rho2];
eos.G    = [G0 G1 G2 G3];
eos.K    = [K0 K1 K2 K3];
eos.a    = [a1 a2 a3];
eos.e    = [e1 e2   ];
eos.rhoc = rhoc;
eos.rhom = rhom;


%{
fprintf('K0    = %16.12e  \n',K0);
fprintf('rhoc  = %16.12e  \n',rhoc);
fprintf('rho0 = %16.12e  \n',rho0);
fprintf('rho1 = %16.12e  \n',rho1);
fprintf('rho2 = %16.12e  \n',rho2);
%}
