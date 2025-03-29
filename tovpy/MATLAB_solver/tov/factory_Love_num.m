% Script to compute Love numbers for various EOS and models


% Set options

% Central density range
rhoc = [5.0e-4 2e-3];
%rhoc = [0.8 6]*1e-3; 
Nmodels = 400;

% l-values and other pars
%lvals  = 2;
lvals  = [2 3 4];

% TOV integration pars
r0 = 1e-12; rN = 25; N = 12000;
rspan = [r0:(rN-r0)/(N-1):rN];
pmin = 0;
verbose = 1;
writefile = 1;  

% choose EOS
k=0;

% k=k+1; eos(k) = PWP4EOS('ALF2');    
% k=k+1; eos(k) = PWP4EOS('BGN1H1');    
% k=k+1; eos(k) = PWP4EOS('ENG');        
% k=k+1; eos(k) = PWP4EOS('FPS');        
% k=k+1; eos(k) = PWP4EOS('H4');
% k=k+1; eos(k) = PWP4EOS('MPA1');
% k=k+1; eos(k) = PWP4EOS('MS1');
% k=k+1; eos(k) = PWP4EOS('MS1b');
% k=k+1; eos(k) = PWP4EOS('2B');
 k=k+1; eos(k) = PWP4EOS('2H');    
% k=k+1; eos(k) = PWP4EOS('HB');
% k=k+1; eos(k) = PWP4EOS('SLy');
    
%{
k=k+1;  eos(k).type = 'poly'; eos(k).name = 'poly'; 
eos(k).K = 123.6489; eos(k).G = 2.;
eos(k).rhoc = +1.229386394647e-03;
%}

% fixme:
%k=k+1;  eos(k).type = 'tab'; eos(k).name = 'APR';   eos(k).file = 'APR.dat';
%k=k+1;  eos(k).type = 'tab'; eos(k).name = 'BSK19'; eos(k).file = 'BSK19.dat';
%k=k+1;  eos(k).type = 'tab'; eos(k).name = 'BSK20'; eos(k).file = 'BSK20.dat';
%k=k+1;  eos(k).type = 'tab'; eos(k).name = 'BSK21'; eos(k).file = 'BSK21.dat';
%k=k+1;  eos(k).type = 'tab'; eos(k).name = 'GNH3';  eos(k).file = 'GNH3.dat';

    
        
% Get down to work    

dr   = diff(rhoc)/(Nmodels-1);
rhoc = rhoc(1):dr:rhoc(end);

for k=1:length(eos)    
    
    fprintf('==> EOS %s\n',eos(k).name)
    
    tov = TOVSequence(rhoc,eos(k),lvals, rspan,pmin, verbose,writefile);    
    
end

