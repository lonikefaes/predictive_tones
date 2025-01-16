function [w, C, F,R] = vb_glm(y,X,W,basis,ard,a0,b0,c0,d0)

% Variational Bayesian optimization routine for general linear model
% of the form y = X w + e, 
%
% Input variables:
%   y         N-by-1 data vector
%   X         N-by-D design matrix 
%   Optional
%   ard       flag for Automatic relevance deteremination (ARD)
%             ard = 1, optimization with ARD (default)
%             ard = 0, optimization without ARD 
%   a0,b0     noninformative inverse-Gamma priors
%   c0,d0     noninformative inverse-Gamma priors
%
% Output variable:
%   w         vector of estimated posterior weights  
%   C         Estimated posterior covariance matrix
%   F         Free energy - lower bound on log-evidence
%
% This code closely follows derivation and implementation by Jan Drugowitsch (2014)
% Variational Bayesian inference for linear and logistic regression
%
% The generative model assumes
%
% p(y | x, w, tau) = N(y | w'x, tau^-1),
%
% with x and y being the rows of the given X and y. w and tau are assigned
% the conjugate normal inverse-gamma prior
%
% p(w, tau | alpha) = N(w | 0, (tau alpha)^-1 I) Gam(tau | a0, b0),
%
% with the hyper-prior
%
% p(alpha) = p(alpha | c0, d0).
%
% The returned posterior parameters (computed by variational Bayesian
% inference) determine a posterior of the form

% N(w1 | w, tau^-1 V) Gam(tau | an, bn).
%
%
% MH2015
% selection of mode. 
if nargin < 3, ard = 1;  end 
% ard = 0 is running VB GLM without ARD,
% ard = 1 is running VB GLM with ARD (default),

%% priors for Gamma distribution (uninformative)
if nargin < 6,  a0 = 1e-2;  end
if nargin < 7,  b0 = 1e-4;  end
if nargin < 8,  c0 = 1e-2;  end
if nargin < 9,  d0 = 1e-4;  end

% pre-process data:

N   = size(W,1);
T   = size(y,2);
W   = kron(basis,W);
X0  = X;
X   = kron(eye(size(basis,2)),X);
[D] = size(X,2);

B   = kron(eye(T),ones(1,N));

y   = y(:);

XtX    = X'*(W'*W)*X;
Xty    = X'*W'*y;

c = c0;
if ard
    d = ones(D,1).*d0;
    p = D;
else
    d = ones(T,1).*d0;
    p = T;
end


% iterate to find hyperparameters
F_new    = -Inf;
Max_iter = 500;

% Loop
for iter = 1:Max_iter
    
    % covariance and weights of linear model
    if ard
        iC  = diag(c./d) + XtX;
    else
        iC  = diag(c./d) + XtX;
    end
    % posterior estimates
    C   =  inv(iC);
    w   =  C*Xty;

    % parameters of noise model
    eet = B*((W*(X*w) - y) .^ 2);

    a   = a0 + N/2; 
    if ard
        b   = b0 + 0.5*(eet + sum(w.^2.*c./d));
    else
        b   = b0 + 0.5*(eet + B*(double(logical(W)*logical(X))*(w.^2)).*c./d);
    end
    % hyperparameters of covariance prior
    if ard
        c   = c0 + 1/2;
        d   = d0 + 0.5 * (a./b.*w.^2 + diag(C));
    else
        c   = c0 + D/2;
        d   = d0 + 0.5 * (a/b*(w'*w) + trace(C));
    end
    
    
    % Free energy without constant terms
    F = - 0.5 * (a./b * eet + sum(sum((W*X).*(W*X*C)))) - 0.5 * logdet(iC) ...
        - b0 * a./b + gammaln(a) - a*log(b) + a ...
        + p*gammaln(c) - c*sum(log(d));
    
    % Free energy must grow (for linear models)!
    if F_new > F
        fprintf('Last bound %6.6f, current bound %6.6f\n', F_new, F);
        error('Free energy should not reduce');
    end
    % stop if change in Free energy is < 0.001%
    if abs(F_new - F) < abs(0.00001 * F)
        fprintf('Converged in %d iterations, F=%6.6f\n', iter, F);
        break
    end
    F_new = F;    
  
end

if ard
    R = diag((d./c)); 
else
    R = (d/c)*eye(D); 
end
if iter == Max_iter
    warning('VB:maxIter','Reached maximum number of iterations.');
end


%% augment Free energy with constant terms

F = F - 0.5 * (N * log(2 * pi) - D) - gammaln(a0) + a0 * log(b0) ...
    + p * (- gammaln(c0) + c0 * log(d0));


