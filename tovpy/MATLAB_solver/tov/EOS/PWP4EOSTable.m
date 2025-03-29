function [nb,e,p] = PWP4EOSTable(eosname,pts)
% PWP4EOSTable.m
% [nb,e,p]=PWP4EOSTable(eosname,N)
% Compute a table of N entries from the PWP4EOS.m output

% constants
mB = 1.675e-24; %1.66*1e-24; %g
cctkMdens_cgs = 6.173895728686583e+17; %6.1764e+17; %g/cm^3
cctkPress_cgs = 5.548820759138184e+38; %5.5511e+38; %dynes/cm^2

% rho
rho_cgs = logspace(0.5,18,pts)';
rho = rho_cgs./cctkMdens_cgs; 

% dimensionless p,e from
% Piecewise Polytropic (4 pieces)
% p = K_i rho^{G_i}
% e = (1+a_i) rho + K_i/(G_i-1) rho^{G_i}

eos = PWP4EOS(eosname); 
rhoi = eos.rho;        
G    = eos.G;
K    = eos.K;
a    = eos.a;
        
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
        
e1 = eos.e(1);
e2 = eos.e(2);        
        
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

% CGS
nb = rho_cgs / mB;
e_cgs = e * cctkMdens_cgs;
p_cgs = p * cctkPress_cgs;


%{ %
figure
subplot 121
loglog(e_cgs,p_cgs,'o-')
xlabel('e'); ylabel('P');
subplot 122
loglog(nb,p_cgs,'o-')
xlabel('n_B'); ylabel('P');
%}

i = [1:pts]';
tab=[i nb e_cgs p_cgs];

%save 'pwp4.dat' tab -ascii

fid = fopen(['pwp_',eosname,'.dat'],'w');
%fprintf(fid,'#\n#\n#\n#\n#\n%d <-- Number of lines\n#\n#\n#\n',pts);
fprintf(fid,'#\n#\n#\n%d\n#\n#\n#\n',pts);
fprintf(fid,'%d %e %e %e\n',tab.');
fclose(fid);

