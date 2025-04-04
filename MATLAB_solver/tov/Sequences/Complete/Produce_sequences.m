GMsun_c2 = 1.476e3; %meters

RHOc = linspace(0.08e-3,9e-3,500);

k = 1;
eos(k).type = 'tab';
eos(k).name = 'BLh';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/BLh_beta_T001.dat';
k = k+1;
eos(k).type = 'tab';
eos(k).name = 'BHBlp';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/BHBlp_beta_T001.dat';
k = k+1;
eos(k).type = 'tab';
eos(k).name = 'DD2';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/DD2_beta_T001.dat';
k = k+1;
eos(k).type = 'tab';
eos(k).name = 'LS200';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/LS220_beta_T001.dat';
k = k+1;
eos(k).type = 'tab';
eos(k).name = 'SFHo';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/SFHo_beta_T001.dat';
k = k+1;
eos(k).type = 'tab';
eos(k).name = 'SLy';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/SLy4_001MeV.dat';
k = k+1;
eos(k).type = 'tab';
eos(k).name = 'TM1';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/TM1_beta_T001.dat';
k = k+1;
eos(k).type = 'tab';
eos(k).name = 'TMA';
eos(k).rhoc = 1.e-3;
eos(k).file = 'EOS/TMA_beta_T001.dat';
m = 1;
eos_pwp(m) = PWP4EOS('H3');
m = m+1;
eos_pwp(m) = PWP4EOS('APR4');
m = m+1;
eos_pwp(m) = PWP4EOS('ALF2');  
m = m+1;
eos_pwp(m) = PWP4EOS('ENG');      
m = m+1;
eos_pwp(m) = PWP4EOS('H4');
m = m+1;
eos_pwp(m) = PWP4EOS('MPA1');
m = m+1;
eos_pwp(m) = PWP4EOS('MS1');
m = m+1;
eos_pwp(m) = PWP4EOS('MS1b');
m = m+1;
eos_pwp(m) = PWP4EOS('2B');
m = m+1;
eos_pwp(m) = PWP4EOS('2H');
m = m+1;
eos_pwp(m) = PWP4EOS('HB');



hold on
for j=1:k
    fprintf('EOS %s\n',eos(j).name)
    tov = TOVSequence(RHOc, eos(j), 2, [], [], 1, 0);

    rhoc = [tov.rhoc];
    M = [tov.M];
    Mb = [tov.Mb];
    R = [tov.R];
    C = [tov.C];
    kl = [tov.kl];
    kl = kl(2:2:end);
    lam = 2./3.*kl.*C.^(-5);
    
    plot(R*GMsun_c2*1e-3, M, 'DisplayName', eos(j).name)
    A = [rhoc; M; Mb; R; C; kl; lam];
    fname = sprintf('Sequences/%s_sequence.txt', eos(j).name);
    fileID = fopen(fname,'w');
    fprintf(fileID, '%15s%15s%15s%15s%15s%15s%15s\r\n','rhoc','M','Mb','R','C','kl','lam');
    fprintf(fileID, '%15.6e%15.6e%15.6e%15.6e%15.4e%15.6e%15.6e\n',A);
    fclose(fileID);
end

for j=1:m
    fprintf('EOS %s\n',eos_pwp(j).name)
    tov = TOVSequence(RHOc, eos_pwp(j), 2, [], [], 1, 0);

    rhoc = [tov.rhoc];
    M = [tov.M];
    Mb = [tov.Mb];
    R = [tov.R];
    C = [tov.C];
    kl = [tov.kl];
    kl = kl(2:2:end);
    lam = 2./3.*kl.*C.^(-5);
    
    plot(R*GMsun_c2*1e-3, M, 'DisplayName',eos_pwp(j).name)
    A = [rhoc; M; Mb; R; C; kl; lam];
    fname = sprintf('Sequences/%s_sequence.txt', eos_pwp(j).name);
    fileID = fopen(fname,'w');
    fprintf(fileID, '%15s%15s%15s%15s%15s%15s%15s\r\n','rhoc','M','Mb','R','C','k2','Lambda2');
    fprintf(fileID, '%15.6e%15.6e%15.6e%15.6e%15.4e%15.6e%15.6e\n',A);
    fclose(fileID);
end


xlabel('$R[km]$','Interpreter','latex')
ylabel('$M[M_\odot]$','Interpreter','latex')
