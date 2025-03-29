function [kappaT, kappaA,kappaB] = TidalCoupConst(MA,CA,kAl,MB,CB,kBl,l)

%TidalCoupConst Compute tidal coupling constant for a binary
%
%   KT = TidalCoupConst(MA,CA,kAl, MB,CBkBl, l)
%   A,B refer to star A and star B
%   M mass 
%   C compactness 
%   k(l) multipolar Love numbers, l=1...l_max 
%   (k(l) and l must be same size!)
%

%{
if (MB<MA)
   tmp = MA;   MA = MB;    MB = tmp;
   tmp = CA;   CA = CB;    CB = tmp;
   tmp = kAl;  kAl = kBl;  kBl = tmp;
   warning('switching A<->B, want MB>MA !')
end
% implementation is symmetric !
%}

% Eq(25)  http://arxiv.org/abs/0911.5041
M         = MA+MB;
MApow2l   = MA.^(2*l);
CAMpow2l1 = (CA*M).^(2*l+1);
MBpow2l   = MB.^(2*l);
CBMpow2l1 = (CB*M).^(2*l+1);

kappaA = 2 * kAl .* MB .* MApow2l ./ CAMpow2l1;
kappaB = 2 * kBl .* MA .* MBpow2l ./ CBMpow2l1;

kappaA(1) = 0; % dummy
kappaB(1) = 0;

kappaT(l) = kappaA + kappaB; 
