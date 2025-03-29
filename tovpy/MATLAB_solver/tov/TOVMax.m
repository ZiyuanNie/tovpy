function tov = TOVMax(rhoc,eos, targets, varargin)

%TOVMax Compute maximum mass of a spherical star sequence
%
%   tov = TOVMax(rhoc,eos, targets)
%   Search for the max of parameter "targets" starting from rhoc  
%   "targets"  is a field in tov struct 
%   rhoc is a guessed central density.
%
%   tov = TOVMax(..., rspan,pmin,Nsurf)
%   Specify options for TOV solver
%
%   tov = TOVMax(..., [],[],pr)
%   Print info
%
%   tov = TOVMAx(..., [],[],[],sv)
%   Save data
%


% Manage args in
rspan = [1e-7:(20-1e-7)/2000:20];
pmin  = 1e-24;
pr    = 1;
sv    = 0;
Nsurf = 4;
if (length(varargin)>5)
    error('too many input args')
end
optargs = {rspan pmin Nsurf pr sv};
newvals = cellfun(@(x) ~isempty(x), varargin);
optargs(newvals) = varargin(newvals);
[rspan,pmin, pr,sv] = optargs{:};


tov = struct([]);


% Check
tov = TOVL(rhoc,eos,2, rspan,pmin,pr,sv,Nsurf);        
if ~isfield(tov,targets)
    error(' %s is not a field in TOV struct',targets);
end


if pr
    fprintf('===> TOV maximum %s\n',targets);
    fprintf(' guess = %.6e\n',tov.(targets));
    tottime = tic;
end


% Compute min
try
    rhoc = fminsearch(@(x) minfun(x, eos,rspan,pmin, targets), rhoc);
    tov = TOVL(rhoc,eos,2, rspan,pmin,0);    
catch ME
    ME
end

if pr
    fprintf(' max   = %.6e\n',tov.(targets));
    fprintf('End %4.3f sec\n',toc(tottime));
end




function x = minfun(rhoc, eos,rspan,pmin, targets)

% Minimum function, assume the quantity targets is positve definite
% (as most of TOV quantities)

tov = TOVL(rhoc,eos,2,rspan,pmin,0);
x   = - abs(tov.(targets));

