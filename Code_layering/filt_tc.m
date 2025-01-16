function [Y] = filt_tc(y,TR,cut_off,center)

if nargin < 4
   center = 0; 
end

[T k] = size(y);

n       = fix(2*(T*TR)/cut_off + 1);
X0      = spm_dctmtx(T,n);
if ~center
    X0  = X0(:,2:end);
end
for i = 1:k
beta    = X0\y(:,i);

Y(:,i)  = y(:,i) - X0*beta;
end