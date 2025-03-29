clear all

ell = [2, 3, 4];
rspan = [1e-7:(20-1e-7)/2000:20];
pmin  = 1e-15;
verbose = 1;
save = 1;

% eos_path = "/home/dur566/Data/BNS/THC/DB/EOS";
% tov_path = "/home/dur566/Data/BNS/THC/DB/EOS";
eos_path = "/home/dur566/Code/TOV/TOVL/EOS/";
tov_path = "/home/dur566/Code/TOV/TOVL/Sequences/Complete/";
eos_name = ["DD2F"];
% eos_name = ["BHBlp_beta_T001", "BLh_beta_T001", "BLQ", ...
%             "DD2_beta_T001", "DD2F_SF1_T01", "LS220_beta_T001", ...
%             "SFHo_beta_T001", "SLy4_001MeV", "DD2F"];

for i = 1:length(eos_name)
    eos.type = 'tab';
    eos.file = eos_path + eos_name(i) + ".dat";
    eos.name = eos_name(i);
    
    N_rhoc = 200;
    rhoc = linspace(1e-4, 5e-3, N_rhoc);
    tov = TOVSequence(rhoc, eos, ell, rspan, pmin, verbose, 0);
    
    if save
        fname = sprintf(tov_path + eos_name(i) + "_sequence.txt");
    
        rc = [tov.rhoc].';
        M  = [tov.M].';        
        Mb = [tov.Mb].';
        R  = [tov.R].';
        C  = [tov.C].';

        n_reshape = length(M);
        kll = reshape([tov.kl],[],n_reshape).';

        T = table(rc,M,Mb,R,C,kll, 'VariableNames', {'rhoc', 'M', 'Mb', 'R', 'C', 'kl'});
        writetable(T,fname, 'Delimiter',' '); 
    end
end