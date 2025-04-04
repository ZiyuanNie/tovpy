% Various examples for computing TOV background, Love numbers, and the
% tidal structure neded for the EOB tidal run
% uncomment what you need


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 1 
    % Compute a star model
    % TOVL.m

    %{  
    % poly EOS
    eos.type = 'poly';
    eos.name = 'poly';
    eos.K    = 100;
    eos.G    = 2;
    ell      = [2 3 4];    
   
    tov = TOVL(1.28e-3,eos,ell,[],[],1);
    %}

    % piecewise-poly EOS
    eos = PWP4EOS('ALF2');    
    %eos = PWP4EOS('BGN1H1');    
    %eos = PWP4EOS('ENG');        
    %eos = PWP4EOS('FPS');        
    %eos = PWP4EOS('H4');
    %eos = PWP4EOS('MPA1');
    %eos = PWP4EOS('MS1');
    %eos = PWP4EOS('MS1b');
    %eos = PWP4EOS('2B');
    %eos = PWP4EOS('2H');    
    %eos = PWP4EOS('HB');
    %eos = PWP4EOS('SLy');
    
    ell = [2 3 4];    
    tov = TOVL(eos.rhoc,eos,ell,[],[],1);
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 2
    % Compute a star with given baryon mass 
    % TOVIterate.m

    Mb = 1.625; % target
    %Mb = 1.44476; % target
    %Mb = 1.77917; % target
    %Mb = 1.90087;
    
    eos.type = 'poly';
    eos.name = 'poly';
    eos.K    = 123.6489; 
    eos.G    = 2;
    eos.rhoc = +1.229386394647e-03;
    ell      = 2;    
   
    % Bisection algorithm    
    tol   = 1e-9;    
    itmax = 200;
    a     = 0.4e-3;      % 1initial guess, bracket solution
    b     = 1.7e-3;            
    tov   = TOVIterate([a b],eos,ell, 'Mb',Mb, tol,itmax);
    
    % Uncontrained zero-finding algorithm (faster)
    tov = TOVIterate(b,eos,ell, 'Mb',Mb, tol,itmax);
    
    fprintf('%.8e',tov.kl(2));
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 3
    % Compute stars with given grav mass and different EOS
    % TOVIterate.m
    
    Mtarget = 1.35;
    
    k=0;
    
    k=k+1; eos(k) = PWP4EOS('ALF2');    
    k=k+1; eos(k) = PWP4EOS('BGN1H1');    
    k=k+1; eos(k) = PWP4EOS('ENG');        
    k=k+1; eos(k) = PWP4EOS('FPS');        
    k=k+1; eos(k) = PWP4EOS('H4');
    k=k+1; eos(k) = PWP4EOS('MPA1');
    k=k+1; eos(k) = PWP4EOS('MS1');
    k=k+1; eos(k) = PWP4EOS('MS1b');
    k=k+1; eos(k) = PWP4EOS('2B');
    k=k+1; eos(k) = PWP4EOS('2H');    
    k=k+1; eos(k) = PWP4EOS('HB');
    k=k+1; eos(k) = PWP4EOS('SLy');
    
    %{ 
    k=k+1;
    eos(k).type = 'poly';
    eos(k).name = 'poly';
    eos(k).K    = 123.6489; 
    eos(k).G    = 2;
    eos(k).rhoc = +1.229386394647e-03;
    %}
    
    for k=1:length(eos)        
        fprintf('==> EOS %s\n',eos(k).name)
        tov = TOVIterate(eos(k).rhoc,eos(k),2, 'M',Mtarget);  
        fprintf(' %.6e\n',tov.kl(2));
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 4
    % Compute a star with maximum mass    
   
    %{ 
    % polytropic EOS
    eos.type = 'poly';
    eos.name = 'poly';
    eos.K    = 100; 
    eos.G    = 2;
    eos.rhoc = 1.28e-3;        

    %eos.type = 'poly'; 
    %eos.name = 'poly';
    %eos.K    = 123.6489; 
    %eos.G    = 2;
    %eos.rhoc = +1.229386394647e-03;
    
    tov = TOVMax(eos.rhoc,eos,'M');
    %}
    
    % piecewise-poly EOS
    k=0;
    k=k+1; eos(k) = PWP4EOS('ALF2');    
    k=k+1; eos(k) = PWP4EOS('BGN1H1');    
    k=k+1; eos(k) = PWP4EOS('ENG');        
    k=k+1; eos(k) = PWP4EOS('FPS');        
    k=k+1; eos(k) = PWP4EOS('H4');
    k=k+1; eos(k) = PWP4EOS('MPA1');
    k=k+1; eos(k) = PWP4EOS('MS1');
    k=k+1; eos(k) = PWP4EOS('MS1b');
    k=k+1; eos(k) = PWP4EOS('2B');
    k=k+1; eos(k) = PWP4EOS('2H');    
    k=k+1; eos(k) = PWP4EOS('HB');
    k=k+1; eos(k) = PWP4EOS('SLy');
    
    for k=1:length(eos)        
        fprintf('==> EOS %s\n',eos(k).name)
        [tov,me] = TOVMax(eos(k).rhoc,eos(k),'M');
        %fprintf(' %.6e %.12e\n',tov.rhoc, log10(tov.rhoc/1.619100425158887e-18));
    end
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 5
    % Compute a sequence use TOVSequence.m

    eos.type = 'poly';
    eos.name = 'poly';
    eos.K    = 100;
    eos.G    = 2;
    ell      = 2;
    
    %rhoc = linspace(0.08e-3,9e-3, 50);            
    rhoc = logspace(-4.0969, -2.0458, 40);            
    tov  = TOVSequence(rhoc,eos,ell, [],[], 1,0);
    
    Mb = [tov.Mb];
    M  = [tov.M];
    R  = [tov.R];    
    kl = reshape([tov.kl],2,40);
    
    figure
    plot(R,M,'o-')    
    xlabel('M'); ylabel('R');
    figure
    plot(rhoc,M./R,'o-')    
    xlabel('\rho_c'); ylabel('C');
    figure
    plot(R,M./R,'o-')    
    xlabel('R'); ylabel('C');
    figure
    plot(M./R,kl(2,:),'o-')    
    xlabel('C'); ylabel('k_2'); 
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 6
    % Compute Love numbers for star sequences
    tic; 
    LoveFactory;
    toc;
        
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 7
    % Prepare tidal struture for EOB run input
    
    Topt.TidalUse = 'yes';
    Topt.PNTidal  = 'nnlo'; 
 
    % object A
    Topt.MA   = 1.51483;
    Topt.RA   = 10.82065;
    Topt.kAl  = [0 0.07890 0 0]; % multipolar tidal parameters (index 1 is dummy)
    Topt.shapehAl  = [0 0.8699 0 0];  % multipolar tidal shape parameters (index 1 is dummy)
     
    % Same for object B
    Topt.MB = 1.51483;
    Topt.RB = 10.82065;
    Topt.kBl = [0 0.07890 0 0];
    Topt.shapehBl  = [0 0.8699 0 0];
            
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
    % Example 8
    % Test tab EOS
    
    ell = [2];    

    eos.rhoc = 1.4002410000000000e-03;
    %{ 
    eos = PWP4EOS('SLy');
    eos.rhoc = 1.4002410000000000e-03
    tovA = TOVL(eos.rhoc,eos,ell,[],[],1);
    %}
    %{ 
    eos.type = 'tab'; 
    eos.name = 'SLy4';  
    eos.file = 'EOS/SLy4.dat';
    tovB = TOVL(eos.rhoc,eos,ell,[1e-8:(10-1e-8)/10000:10],1e-17,1);
    %}

    %{ % 
    eos.type = 'poly';
    eos.name = 'poly';
    eos.K    = 100;
    eos.G    = 2;
    tovA = TOVL(1.28e-3,eos,ell,[],[],1);
    %}

    eos.type = 'tab'; 
    eos.name = 'polyTAB';  
    eos.file = 'EOS/polytropic.d';
    tovB = TOVL(1.28e-3,eos,ell,[],[],1);

end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
    
    % Example 9
    % Compute M vs R for a sequences of TOV w\ tab EOS

    eos.type = 'tab';
    ell      = 2;
    
    rhoc = linspace(0.08e-3,9e-3, 50);            
    %rhoc = logspace(-4.0969, -2.0458, 40);   

    eos.type = 'tab'; 
    
    eos_file = {'SLy4_001MeV.dat', 'DD2_beta_T001.dat','BHBlp_beta_T001.dat','LS220_beta_T001.dat','SFHo_beta_T001.dat'};
    eos_name = {'SLy4','DD2','BHBlp','LS220','SFHo'};  
    rhoc  = { linspace(5e-4,6e-3,50), ...
        linspace(6e-4,2.75e-3,50), ...
        linspace(6e-4,2.75e-3,50), ...
        linspace(5e-4,5e-3,50), ...
         linspace(5e-4,5e-3,50), ...
        };
    
    figure
    hold all
    
    for i = 1:length(eos_name)

        eos.file = ['EOS/',eos_file{i}];
        eos.name = eos_name{i};
        %tov  = TOVSequence(rhoc{i},eos,ell, [],1e-24, 1,0);
        tov  = TOVSequence(rhoc{i},eos,ell, [],[], 1,0);
  
        M  = [tov.M];
        R  = [tov.R];    
        plot(R,M,'-')    
  
    end
    legend(eos_name)
    ylabel('M'); xlabel('R');
    
end
