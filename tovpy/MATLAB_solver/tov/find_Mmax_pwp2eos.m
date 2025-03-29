% Script to find maximum mass and M=1.35 stars of 2-segments PWP EOS
% Compare with Lackey et al. 1109.3402 Tab. I 

k=0;
k=k+1; eos(k) = PWP2EOS('p.3G2.4_Bss');
% k=k+1; eos(k) = PWP2EOS('p.3G2.7_Bs');
% k=k+1; eos(k) = PWP2EOS('p.3G3.0_B');
% k=k+1; eos(k) = PWP2EOS('p.3G3.3');
% k=k+1; eos(k) = PWP2EOS('p.4G2.4_HBss');
% k=k+1; eos(k) = PWP2EOS('p.4G2.7_HBs');
% k=k+1; eos(k) = PWP2EOS('p.4G3.0_HB');
% k=k+1; eos(k) = PWP2EOS('p.4G3.3');
% k=k+1; eos(k) = PWP2EOS('p.5G2.4');
% k=k+1; eos(k) = PWP2EOS('p.5G2.7');clear al
% k=k+1; eos(k) = PWP2EOS('p.5G3.0_H');
% k=k+1; eos(k) = PWP2EOS('p.5G3.3');
%k=k+1; eos(k) = PWP2EOS('p.6G2.4'); % <---- fixme
% k=k+1; eos(k) = PWP2EOS('p.6G2.7');
% k=k+1; eos(k) = PWP2EOS('p.6G3.0');
% k=k+1; eos(k) = PWP2EOS('p.6G3.3');
%k=k+1; eos(k) = PWP2EOS('p.7G2.4'); % <---- fixme
% k=k+1; eos(k) = PWP2EOS('p.7G2.7');
% k=k+1; eos(k) = PWP2EOS('p.7G3.0_1.5H');
% k=k+1; eos(k) = PWP2EOS('p.7G3.3');
% k=k+1; eos(k) = PWP2EOS('p.9G3.0_2H');


for k=1:length(eos)
    fprintf('==> EOS %s\n',eos(k).name)
    tov = TOVMax(eos(k).rhoc,eos(k),'M');
    %tov = TOVMax(eos(k).rhoc,eos(k),'M',[1e-9:14/2000:14],1e-22);
    fprintf('rhom = %.6e;\nlgrhom = %.12e;\n',tov.rhoc, log10(tov.rhoc/1.619100425158887e-18));
end


for k=1:length(eos)
    fprintf('==> EOS %s\n',eos(k).name)
    tov = TOVIterate(eos(k).rhoc,eos(k),2, 'M', 1.35);  
    fprintf('lgrhoc = %.12e; ',log10(tov.rhoc/1.619100425158887e-18));
    fprintf('R[km] = %.2f;\n',tov.R*1.476701332464468);
end
