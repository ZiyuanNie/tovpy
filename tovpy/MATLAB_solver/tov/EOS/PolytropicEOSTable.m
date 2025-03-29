function [nb,e,p] = PolytropicEOSTable(k,g,pts)
% PolytropicEOSTable.m
% [nb,e,p]=PolytropicEOSTable(k,g,N)
% Compute a table of N entries from the EoS

% constants
mB = 1.675e-24; %1.66*1e-24; %g
cctkMdens_cgs = 6.173895728686583e+17; %6.1764e+17; %g/cm^3
cctkPress_cgs = 5.548820759138184e+38; %5.5511e+38; %dynes/cm^2

% rho
rho_cgs = logspace(0.5,18,pts)';
rho = rho_cgs./cctkMdens_cgs; 

% dimensionless p,e
p = k*rho.^g;
e = rho + p./(g-1);

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

%save 'polytropic.dat' tab -ascii

fid = fopen('polytropic.d','w');
%fprintf(fid,'#\n#\n#\n#\n#\n%d <-- Number of lines\n#\n#\n#\n',pts);
fprintf(fid,'#\n#\n#\n%d\n#\n#\n#\n',pts);
fprintf(fid,'%d %e %e %e\n',tab.');
fclose(fid);

