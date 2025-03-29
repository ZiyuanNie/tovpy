% Script to convert from Lorene to TOV.m format [fm->cm]

% DD2
%eosfile = {'Hempel_DD2_eostable_NFBetaEq_pnu_T=00.500.lorene'};
%outfile = {'DD2_beta_T005.dat'};

%eosfile = {'DD2-0.01MeV_22-Jun-2018.lorene'};
%outfile = {'DD2_beta_T001.dat'};

% LS220
%eosfile = {'LS220B0.lorene'};
%outfile =  {'LS220_beta_T001.dat'};

% SFHo
%eosfile = {'SFHo.lorene'};
%outfile = {'SFHo_beta_T005.dat'};

% BHBlp
%eosfile = {'BHBlp.dat'}; % wrong
eosfile = {'bhb_lp_25-sept-2017.lorene'};
outfile = {'BHBlp_beta_T005.dat'};

%eosfile = {'BSK19_newP.dat' 'BSK20_newP.dat' 'BSK21_newP.dat' 'GNH3.dat' 'SLy4.dat'}
%outfile = {'BSK19_newP.dat' 'BSK20_newP.dat' 'BSK21_newP.dat' 'GNH3.dat' 'SLy4.dat'}
%eosfile = {'BSK19.dat' 'BSK20.dat' 'BSK21.dat'}
%outfile = {'BSK19.dat' 'BSK20.dat' 'BSK21.dat'}


for k=1:length(eosfile)

% go
aux = LoadEOSTab( eosfile{k} );
N = length(aux{2});
out_table = zeros(N,4);
out_table(:,1) = 1:N;
out_table(:,2) = aux{2}*1e39;
out_table(:,3) = aux{3};
out_table(:,4) = aux{4};

fid = fopen(outfile{k},'w');
fprintf(fid,'#\n#\n#\n%d <-- Number of lines\n#\n#\n#\n',N);
fprintf(fid,'%d %e %e %e\n',out_table.');
fclose(fid);


end



function col = LoadEOSTab(fname)
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
end