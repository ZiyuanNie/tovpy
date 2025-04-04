function tov = TOVSequence(rhoc,eos,ell, varargin)

%TOVSequence Compute sequence of spherical stars and Love numbers
%
%   tov = TOVSequence(rhoc,eos,ell)
%   Compute a TOV sequence for central densities in rhoc 
%
%   tov = TOVSequence(..., rspan,pmin)
%   Specify options for TOV solver
%
%   tov = TOVSequence(..., [],[],pr)
%   Print info
%
%   tov = TOVSequence(..., [],[],[],sv)
%   Save data
%


% Manage args in
rspan = [1e-7:(20-1e-7)/2000:20];
pmin  = 1e-24;
pr    = 1;
sv    = 0;
if (length(varargin)>4)
    error('too many input args')
end
optargs = {rspan pmin pr sv};
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[rspan,pmin, pr,sv] = optargs{:};


if pr
    fprintf('===> TOV sequence\n');
    tottime = tic;
end

% Compute sequence
computeLove = sum((ell>1));
c = 0;
for k=1:length(rhoc)
    
    try
        
        aux = TOVL(rhoc(k),eos,ell,rspan,pmin, 0);    
        
        if computeLove & ...
                (sum(isnan(aux.kl)) || sum(isnan(aux.hl)) || sum(isnan(aux.jl)))
            % some calculation may fail: exclude them
            continue;
        end
        
        c = c+1;
        tov(c) = aux;
        
        if pr
            fprintf('====> rhoc(%d) = %.4e',k,rhoc(c));
            fprintf(' M = %.4e C = %.6e', tov(c).M,tov(c).C);
            if computeLove
                fprintf(' k_2 = %.6e',tov(c).kl(2));
            end
            fprintf('\n');
        end
    
    catch ME
        ME
    end 
    
end

if sv
    % Save data files
    fname = sprintf('TOVSeq_EOS%s.dat',eos.name);

    rc = [tov.rhoc].';
    Mb = [tov.Mb].';
    M  = [tov.M].';
    R  = [tov.R].';
    Rp = [tov.Rp].';
    C  = [tov.C].';   
         
    if computeLove

        kl = reshape(real([tov.kl]),length(ell)+1,length(tov)); 
        hl = reshape(real([tov.hl]),length(ell)+1,length(tov));    
        jl = reshape(real([tov.jl]),length(ell)+1,length(tov));    
    
        kl = kl(2:end,:).';
        hl = hl(2:end,:).';
        jl = jl(2:end,:).';

        headl = sprintf('rhoc C Mb M R Rp kl hl jl (l=%s)',mat2str(ell));
        WriteASCII(fname,rc,[C Mb M R Rp kl hl jl], headl);
    
    else 

        headl = sprintf('rhoc C Mb M R Rp C');
        WriteASCII(fname,rc,[C Mb M R Rp], headl);
   
    end

end


if pr
    fprintf('End %4.3f sec\n',toc(tottime));
end
