function pval = rayleigh_pval(r, n);
% rayleigh_pval - confidence interval p of Rayleigh test
% pval = rayleigh_pval(r, N)
%   r is vector strength
%   N is # events
%   pval is p value.
%
%   Inputs r N may be arrays having compatible sizes (see SameSize).
%   Code adapted from Philipp Berens toolbox on circular stats.
%
%   See also rayleigh_test, rayleigh_rmin.

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

[r, n] = SameSize(r, n);

R = n.*r;
% compute Rayleigh's z (equ. 27.2)
z = R.^2 ./ n;
% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n.^2-R.^2))-(1+2*n));

