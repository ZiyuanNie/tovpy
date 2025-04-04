function tov = TOVL(rhoc, eos, ell, varargin)

%TOVL Solve TOV for spherical equilibrium star configuration
%   and perturbation equation to compute Love numbers
%
%   Schwarzschild coordinates:
%     ds^2 = -exp(2 nu)dt^2 + exp(2 lambda)dr^2 + r^2 dOmega^2(theta,phi)
%   Units G=c=Msun=1 
%
%   tov = TOVL(rhoc, eos, ell) 
%
%   Given
%    rhoc   : central rest-mass density
%    eos    : a structure to give the EOS (see below)
%    ell    : multipolar index(es) of the pertrubation
%
%   Return a data structure with 
%    M   : gravitational mass
%    Mb  : baryonic mass
%    R   : radius
%    Rp  : proper radius
%    C   : compactness
%    r   : radial grid
%    m   : mass function, m(r) = 0.5 r (1 - 1/g_11)
%    e   : total energy density
%    rho : rest-mass density
%    p   : pressure
%    cs2 : square of sound speed
%    nu  : metric function, nu(r) = 0.5 ln(-g_00)
%    lam : metric function, lambda(r) = -0.5 ln(g_11) 
%
%    H   : even parity pertrubation
%    Psi : odd    "        "
%    yr  : = r H'(r)/H(r)
%    yeven  : R H'(R)/H(R)
%    yodd   : R Psi'(R)/Psi(R)
%    kl  : Love number (apsidal constant), size(ell)
%    hl  : shape Love number 
%    jl  : Love number
%
%  EOS structure has 3 possibilities
%
%   - Polytropic EOS
%   eos.name = 'poly' 
%   eos.K = polytropic constant
%   eos.G = polytropic index
%
%   - Piecewise-Polytropic EOS (4 segments)
%   eos.name = 'pwp4'
%   eos.K0   = low-density polytropic constant
%   eos.Gi   = polytropic indexes i=1...4
%   eos.rhoi = rho-interfaces i=1...4
%
%   - Tabulated EOS
%   eos.name = 'tab'
%   eos.file = filename
%   The format required is ASCII 3 cols:
%       1: number density n [1/cm^3]
%       2: energy density e [g/cm^3] 
%       3: pressure       p [dyn/cm^2] 
%
%
%   tov = TOVL(..., rspan ) 
%   Specify radial grid {[rmin rmax],[rmin:dr:rmax]}
%
%   tov = TOVL(...,[], pmin ) 
%   Specify minimum pressure (1e-24)
%
%   tov = TOVL(...,[],[], pr ) 
%   Screen info
%
%   tov = TOVL(...,[],[],[], sv ) 
%   Save matab file
%
%   tov = TOVL(...,[],[],[],[],Nsurf ) 
%   number of additional Euler steps for star surface computation
%
%   Refs
%   C.V.Misner, K.S.Thorne & J.A.Wheeler, "Gravitation", Freeman, 1974, Chap.23
%   E.N.Glass & L.Lindblom, ApJS 53:93-103, 1983
%   T.Hinderer, Astrophys.J.677:1216-1220, 2008 http://arxiv.org/abs/0711.2420
%



% Manage args in
rspan = [1e-9:(20-1e-9)/2000:20];
pmin  = 0; 
pr    = 1;
sv    = 0;
Nsurf = 4;
if (length(varargin)>5)
    error('too many input args')
end
optargs = {rspan pmin pr sv Nsurf};
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[rspan, pmin, pr, sv, Nsurf] = optargs{:};


if pr
    fprintf('===> TOV run\n');
    tottime = tic;
end


% Set options
odeoptions = odeset('RelTol',1e-10,'AbsTol',1e-12,...
    'Refine',6,...
    'Events',@(r,x) PresZero(r,x,pmin));

if strcmp(eos.type,'tab')
    tab = LoadEOSTab(eos.file);
    eos.tab = tab;
    clear tab;
end

dr = min(diff(rspan));

% Intial data
[p0,e0] = EOS_r2p(rhoc, eos.type, eos);
rmin    = rspan(1);
m0      = 4/3*pi*e0*rmin^3;

% Integrate
tov.ell=ell;
for l=ell
    
    % Initial data for l-perturbation
    a0    = 1;
    H0    = a0 * rmin^l;
    dH0   = a0 * l*rmin^(l-1);
    Psi0  = a0 * rmin^(l+1);
    dPsi0 = a0 * (l+1)*rmin^l;
    X0    = [m0 p0 0 H0 dH0 Psi0 dPsi0];

    % Integration
    if pr
        fprintf('====> EOS = %s rhoc = %.6e l = %d ...\n',eos.name,rhoc,l);
        tstart = tic;           
    end
    [r,X,rstar,XE,IE] = ode113(@(r,x) TOVODE(r,x, eos,l), ...
        rspan,X0,odeoptions);
    if pr
        fprintf(' Done %4.3f sec\n',toc(tstart));
    end

    % Checks
    if isempty(rstar)
        error('Star radius not reached, rerun with different rspan/pmin.')
    end
    
    % Determine the radius 
    % Zero of pressure is now located at 
    %    X(end-1,2) < P(R) < X(end,2)
    % -------------------------------------------
    % 1. do Nsurf steps with Euler and get to P<0
    %dr = 6e-3;
    XE = X(end-1,:).'; %transpose
    dr = min(diff(rspan));
    for i=1:Nsurf
        %fprintf('%.16e \n',XE)
        XE = XE + dr * TOVODE(r(end-1),XE, eos, l);
        if XE(2)<0.
            break
        end
    end
    if XE(2)>0.
            error('Can''t locate star radius (p=%.16e), rerun with different Nsurf (%d) or change rspan/pmin. ',XE(2),Nsurf)
    end
    % 2. interpolate to P(R)=0 
    rstar  = interp1([X(end-1,2),XE(2)],r(end-1:end),0,'linear');
    % 3. advance to rstar with Euler step 
    dr = rstar - r(end-1);
    X(end,:) = X(end-1,:).' + dr * TOVODE(r(end-1),X(end-1,:).', eos, l);
    raux = r; raux(end) = rstar;
    % 4. regrid    
    r = [raux(1):(raux(end)-raux(1))/(length(raux)-1):rstar].';
    X = interp1(raux,X, r,'linear');
    % Compute Love numbers
    %R    = rstar(1);
    %M    = XE(1,1);
    R    = rstar;
    M    = X(end,1);
    C    = M/R;
    H    = X(:,4);
    dH   = X(:,5);
    Psi  = X(:,6);
    dPsi = X(:,7);
    
    yr          =  r.*dH./H;
    yeve        =  R.*dH(end)./H(end);
    yodd        =  R.*dPsi(end)./Psi(end);
    
    [kell, hell] =  ComputeLove(C,yeve,l, 'even');
    jell        =  ComputeLove(C,yodd,l, 'odd' );

    tov.yr{l}   = yr;
    tov.yeven{l}= yeve;
    tov.yodd{l} = yodd;
    tov.H{l}    = H;
    tov.dH{l}   = dH;
    tov.Psi{l}  = Psi;
    tov.dPsi{l} = dPsi;
    
    kl(l) = kell;
    hl(l) = hell;
    jl(l) = jell;

    if pr
        fprintf('  k(%d)  = %+.12e\n',l,kell);
        fprintf('  h(%d)  = %+.12e\n',l,hell);
        fprintf('  s(%d)  = %+.12e\n',l,jell);
    
    end
      
end

m      = X(:,1); 
tov.kl = kl; 
tov.hl = hl; 
tov.jl = jl; 

% Store background
tov.r    = r;
tov.rhoc = rhoc;
tov.M    = M;
tov.R    = R;
tov.C    = C;
tov.m    = m; 
tov.p    = X(:,2);
tov.nu   = X(:,3) + 0.5*log(1-(2*M)/R)-X(end,3); % Schwarzschild match 
tov.lam  = - 0.5*log(1-(2*m)./r); 

[e,rho,cs2] = EOS_p2e(X(:,2), eos.type,eos);
tov.e   = e;
tov.rho = rho;
tov.cs2 = cs2;


% Baryonic mass
tov.Mb = trapz( r, 4*pi*r.^2.*rho./sqrt(1-2*m./r) );

% Proper radius
tov.Rp = trapz( r, 1./sqrt((1-2*m./r)) );


if sv
    % filename ( eos_mass )
    filename = strcat(eos.name,sprintf('_M%4.3f',M));
    save(filename,'tov','eos');
end

if pr
     PrintSFields(tov,1);    
     fprintf('End %4.3f sec\n',toc(tottime));
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dx = TOVODE(r,x, eos, ell)     
% Tolman-Oppenheimer-Volkoff eqs + perturbation eqs

dx = zeros(7,1);       

% TOV eqs
m   = x(1); % mass included in a sphere of radius r
p   = x(2); % pressure
r2  = r^2;
r3  = r^3;
tmp = (m + 4*pi*r3*p)./(r*(r-2*m));

[e,rho,cs2] = EOS_p2e(p, eos.type,eos);        

dx(1) = 4*pi*e*r2;                                   
dx(2) = -(e+p)*tmp;
dx(3) = tmp;             

if ell>1
    
    % Perturbation eqs
    exp_lam = 1./(1-2*m/r);
    dnudr2  = (2*tmp)^2;

    ceve0 = - ( 2/r + exp_lam*(2*m/r2 + 4*pi*r*(p-e)) );
    ceve1 = - exp_lam*(-ell*(ell+1)/r2 + 4*pi*(5*e + 9*p + (e+p)/cs2)) + dnudr2;

    codd0 =  exp_lam*( ell*(ell+1)/r2 - 6*m/r3 + 4*pi*(e-p) );  
    codd1 = -exp_lam*( 2*m+4*pi*r3*(p-e) )./r2;

    dx(4) = x(5);
    dx(5) = ceve0*x(5) + ceve1*x(4);    
    dx(6) = x(7);
    dx(7) = codd0*x(6) + codd1*x(7);

end


% Event function
function [value, isterminal, direction] = PresZero(r,x, pmin)
%fprintf('%e\n',x(2))
value      = x(2) - pmin;
isterminal = 1;   
direction  = 0; %-1;   


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOS routines
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [e,rho,cs2] = EOS_p2e(p, type, opt)
% Compute e(p), rho(p) cs2(p)

switch type
    
    case 'poly'

        % EOS: Polytropic
        % p = k rho^G

        K = opt.K;
        G = opt.G;
        
        e   = (p>0).*(p/K).^(1/G) + p/(G-1);
        rho = (p/K).^(1/G);                        
        cs2 = K*G*(G-1)./(G-1 + K*G*rho.^(G-1)) .*rho.^(G-1);
    
    case 'pwp4'
        
        % EOS: Piecewise Polytropic (4 pieces)
        % p = K_i rho^{G_i}
        % e = (1+a_i) rho + K_i/(G_i-1) rho^{G_i}
            
        rhoi = opt.rho;        
        G    = opt.G;
        K    = opt.K;
        a    = opt.a;
        
        G0 = G(1);
        G1 = G(2);
        G2 = G(3);
        G3 = G(4);
        
        K0 = K(1);
        K1 = K(2);
        K2 = K(3);
        K3 = K(4);        

        a1 = a(1);
        a2 = a(2);
        a3 = a(3);
        
        e1 = opt.e(1);
        e2 = opt.e(2);        
        
        rho0 = rhoi(1);
        rho1 = rhoi(2);
        rho2 = rhoi(3);
        
        p0 = K0*rho0^G0;
        p1 = K1*rho1^G1;
        p2 = K2*rho2^G2;   

        % init =0
        rho = 0*p;
        e   = rho;
        cs2 = rho;

        I = find(p>p2);
        rho(I) = (p(I)/K3).^(1/G3);
        e(I)   = (1+a3)*rho(I) + K3/(G3-1).*rho(I).^G3;
        cs2(I) = G3.*p(I)./(e(I) + p(I));

        I = find( (p>p1) & (p<=p2) );
        rho(I) = (p(I)/K2).^(1/G2);
        e(I)   = (1+a2)*rho(I) + K2/(G2-1).*rho(I).^G2;
        cs2(I) = G2.*p(I)./(e(I) + p(I));

        I = find( (p>p0) & (p<=p1) );
        rho(I) = (p(I)/K1).^(1/G1);
        e(I)   = (1+a1)*rho(I) + K1/(G1-1).*rho(I).^G1;
        cs2(I) = G1.*p(I)./(e(I) + p(I));

        I = find(p<=p0);
        rho(I) = (p(I)/K0).^(1/G0);
        e(I)   = rho(I) + K0/(G0-1).*rho(I).^G0;
        cs2(I) = G0.*p(I)./(e(I) + p(I));

        I = find(p<=0);
        p(I)   = 0;
        rho(I) = 0;
        e(I)   = 0;
        cs2(I) = 0;

    case 'tab'
        
        % EOS: Table
        lgp = log10(p);
  
        % linear interp
        %TABlgr = opt.tab.lgr;
        %TABlge = opt.tab.lge;
        %TABlgp = opt.tab.lgp;
        %lgr        = lininterp(TABlgp,TABlgr, lgp);
        %[lge,dlge] = lininterp(TABlgp,TABlge, lgp);
  
        % spline interp
        lgr = ppval(opt.tab.PP_pr, lgp);
        lge = ppval(opt.tab.PP_pe, lgp);
        dlge = ppval(opt.tab.PP_peD, lgp);

        e   = 10.^lge;
        rho = 10.^lgr;        
        cs2 = p./e./dlge; % !

    otherwise        
        
        error('unknown EOS type %s',type);        
        
end


function [p,e] = EOS_r2p(rho, type, opt)
% Compute p(rho) and e(rho)

switch type
        
    case 'poly'
        
        % EOS: Polytropic
        % p = k rho^G

        K = opt.K;        
        G = opt.G;
        
        p = K*(rho).^(G); 
        e = rho + p/(G-1);
        
    case 'pwp4'
        
        % EOS: Piecewise Polytropic (4 pieces)
        % p = K_i rho^{G_i}
        % e = (1+a_i) rho + K_i/(G_i-1) rho^{G_i}

        rhoi = opt.rho;        
        G    = opt.G;
        K    = opt.K;
        a    = opt.a;
        
        G0 = G(1);
        G1 = G(2);
        G2 = G(3);
        G3 = G(4);
        
        K0 = K(1);
        K1 = K(2);
        K2 = K(3);
        K3 = K(4);        

        a1 = a(1);
        a2 = a(2);
        a3 = a(3);
        
        e1 = opt.e(1);
        e2 = opt.e(2);        
        
        rho0 = rhoi(1);
        rho1 = rhoi(2);
        rho2 = rhoi(3);           

        % init =0
        p   = 0*rho;
        e   = p;

        I = find(rho>rho2);
        p(I) = K3*rho(I).^G3;
        e(I) = (1+a3)*rho(I) + K3/(G3-1).*rho(I).^G3;
 
        I = find( (rho>rho1 & rho<=rho2) );
        p(I) = K2*rho(I).^G2;
        e(I) = (1+a2)*rho(I) + K2/(G2-1).*rho(I).^G2;
 
        I = find( (rho>rho0 & rho<=rho1) );
        p(I) = K1*rho(I).^G1;
        e(I) = (1+a1)*rho(I) + K1/(G1-1).*rho(I).^G1;

        I = find( (rho<=rho0) );
        p(I) = K0*rho(I).^G0;
        e(I) = rho(I) + K0/(G0-1).*rho(I).^G0;
     
        I = find(rho<=0);
        p(I)   = 0;
        e(I)   = 0;
        
        
    case 'tab'
        
        % EOS: Table        
        lgrho = log10(rho);
        % linear interp
        %TABlgr = opt.tab.lgr;
        %TABlgp = opt.tab.lgp;
        %TABlge = opt.tab.lge;
        %lgp = lininterp(TABlgr,TABlgp, lgrho);
        %lge = lininterp(TABlgr,TABlge, lgrho);
        
        % spline
        lgp = ppval(opt.tab.PP_rp, lgrho);
        lge = ppval(opt.tab.PP_re, lgrho);
        
        p   = 10.^lgp;        
        e   = 10.^lge;
        
    otherwise
        
        error('unknown EOS type %s',type);
        
end
        

% Some routines for EOS tables


function tab = LoadEOSTab(fname)
% Load EOS table
% Format (ASCII, 4 cols):
%
% #
% #
% #
% nlines
% #
% #
% #
% index n[1/cm^3] e[g/cm^3] p[dyn/cm^2]
%

fid = fopen(fname,'r');
col = textscan(fid,'%n%n%n%n','Headerlines',7);
fclose(fid);
  
%Length_fm = 1.476701332464468e+18; % cm

units.vol  = 3.2202e+15;             
units.pres = 5.548820759138184e+38;  
units.mden = 6.173895728686583e+17;  
units.mB   = 1.675e-24;      
units.Msun = 1.9889e+33;

tab.rho = col{2} * units.mB / units.mden; 
%tab.lgn = log10( col{2} * units.vol );
tab.lgr = log10( tab.rho ); 
%tab.lgr = log10( col{2} * units.mB * units.vol / units.mden ); % rho = mB n % !
tab.lge = log10( col{3} / units.mden );
tab.lgp = log10( col{4} / units.pres );

if isempty(tab.lgr)
    error('problem loading EOS table');
end

% for spline store the piecewise interpolant poly (do it once!)
tab.PP_pr = spline(tab.lgp,tab.lgr);
tab.PP_pe = spline(tab.lgp,tab.lge);
T = tab.PP_pe.coefs; % temp storage
T = [3*T(:,1),2*T(:,2),T(:,3)];
tab.PP_peD = mkpp(tab.PP_pe.breaks,T); % derivative !
tab.PP_rp = spline(tab.lgr,tab.lgp);
tab.PP_re = spline(tab.lgr,tab.lge);
