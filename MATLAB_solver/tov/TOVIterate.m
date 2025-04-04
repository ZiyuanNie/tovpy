function tov = TOVIterate(rhoc,eos,ell, targets,targetv, varargin)

%TOVIterate Compute spherical stars configuration with given parameter
%
%   tov = TOVIterate(rhoc,eos,ell, targets,targetv)
%   Compute model with fix value of parameter "target"
%   "targets"  is a field in tov struct 
%   "targetv"   is the fix value
%   the iteration stops if
%   | TOV.(targets) - targetv |/targetv < tol
%   rhoc can be 
%      . an interval of central density that brakets the root 
%        => use bisection algorithm
%      . a guess value 
%        => use unconstrained zero-finding algorithm
%
%   tov = TOVIterate(..., tol,itmax)
%   Specify relative tolerance and max iteration
%
%   tov = TOVIterate(..., [],[], rspan,pmin)
%   Specify options for TOV solver
%
%   tov = TOVIterate(..., [],[], [],[],pr)
%   Print info
%



% Manage args in
tol = 1e-6;
itmax = 200;
rtov = [1e-6:(20-1e-6)/2000:20];
ptovmin = 0;
pr = 0;
if (length(varargin)>5)
    error('too many input args')
end
optargs = {tol itmax rtov ptovmin pr};
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);

[tol,itmax, rtov,ptovmin, pr] = optargs{:};


% Iterate
if pr
    fprintf('===> TOV iterate\n');
    fprintf(' tol   = %.6e\n',tol);
    fprintf(' itmax = %.6e\n',itmax);
    tottime = tic;
end


tov = struct([]);


if (sum(length(rhoc))==2)
    try
        % use bisection
        tov = bisect(rhoc,eos,ell, targets,targetv,tol,itmax, rtov,ptovmin);
    catch ME
        ME
    end    
else
    try
        % use unconstrained minimization
        rhoc = fzero(@(x) fz(x, eos,ell, targets,targetv, rtov,ptovmin), rhoc);
        tov  = TOVL(rhoc,eos,ell, rtov,ptovmin,0);    
    catch ME
        ME
    end    
end


if pr
    PrintSFields(tov,1);
    fprintf('End %4.3f sec\n',toc(tottime));
end




function tc = bisect(rhoab,eos,ell, targets,targetv,tol,itmax, rtov,ptovmin)

% Bisection routine to iterate on parameter "target"
% | TOV.(targets) - targetv |/targetv < tol

% initial guess, bracket solution
a  = rhoab(1);
b  = rhoab(2);

% bisection
ta = TOVL(a,eos,ell, rtov,ptovmin,0);
tb = TOVL(b,eos,ell, rtov,ptovmin,0);
if ~isfield(ta,targets)
    error(' %s is not a field in TOV struct',targets);
end

% check bracketing
fa = (targetv - ta.(targets));
fb = (targetv - tb.(targets));
if sign(fa)==sign(fb)
    error('Root not in bracket');
end

% bisection loop
it  = 0;
err = 1;
while err>tol
    it = it+1;
    if it>itmax
        break;
    end
    c  = a + 0.5*(b-a);
    tc = TOVL(c,eos,ell, rtov,ptovmin,0);
    fc = (targetv - tc.(targets));
    if sign(fa)~=sign(fc); % root lies in [a,c]
        b  = c;
        fb = fc;
    else % root lies in [c,b]
        a  = c;
        fa = fc;
    end
    err = abs(fc)/targetv;
end
if it>itmax || err>tol
    error('max iteration reached (%d,%.6e)',it,tol);
    %warning('max iteration reached (%d,%.6e)',it,tol);
end




function m = fz(rho,eos,ell, targets,targetv, rtov,ptovmin)

% Function for unconstrained nonlinear minimization routine 

t = TOVL(rho,eos,ell, rtov,ptovmin,0); 
m = t.(targets)-targetv;



