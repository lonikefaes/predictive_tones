function [y X dFx] = spm_int_IT(P,M,U)
% integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D;
% FORMAT [y] = spm_int(P,M,U)
% P   - model parameters
% M   - model structure
%   M.delays - sampling delays (s); a vector with a delay for each output
%
% U   - input structure or matrix
%
% y   - response y = g(x,u,P)
%__________________________________________________________________________
% Integrates the bilinear approximation to the MIMO system described by
%
%    dx/dt = f(x,u,P) = A*x + u*B*x + C*u + D
%    y     = g(x,u,P) = L*x;
%
% at v = M.ns is the number of samples [default v = size(U.u,1)]
%
% spm_int will also handle static observation models by evaluating
% g(x,u,P).  It will also handle timing delays if specified in M.delays
%
%--------------------------------------------------------------------------
%
% SPM solvers or integrators
%
% spm_int_ode:  uses ode45 (or ode113) which are one and multi-step solvers
% respectively.  They can be used for any ODEs, where the Jacobian is
% unknown or difficult to compute; however, they may be slow.
%
% spm_int_J: uses an explicit Jacobian-based update scheme that preserves
% nonlinearities in the ODE: dx = (expm(dt*J) - I)*inv(J)*f.  If the
% equations of motion return J = df/dx, it will be used; otherwise it is
% evaluated numerically, using spm_diff at each time point.  This scheme is
% infallible but potentially slow, if the Jacobian is not available (calls
% spm_dx).
%
% spm_int_E: As for spm_int_J but uses the eigensystem of J(x(0)) to eschew
% matrix exponentials and inversion during the integration. It is probably
% the best compromise, if the Jacobian is not available explicitly.
%
% spm_int_B: As for spm_int_J but uses a first-order approximation to J
% based on J(x(t)) = J(x(0)) + dJdx*x(t).
%
% spm_int_L: As for spm_int_B but uses J(x(0)).
%
% spm_int_U: like spm_int_J but only evaluates J when the input changes.
% This can be useful if input changes are sparse (e.g., boxcar functions).
% It is used primarily for integrating EEG models
%
% spm_int: Fast integrator that uses a bilinear approximation to the
% Jacobian evaluated using spm_bireduce. This routine will also allow for
% sparse sampling of the solution and delays in observing outputs. It is
% used primarily for integrating fMRI models (see also spm_int_D)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_int.m 6856 2016-08-10 17:55:05Z karl $
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U), u.u = U; U = u; end
try, dt = U.dt; catch, U.dt = 1; end
 
% number of times to sample (v) and number of microtime bins (u)
%--------------------------------------------------------------------------
u       = size(U.u,1);
try,  v = M.ns;  catch, v = u;   end


% get expansion point
%--------------------------------------------------------------------------
x = [spm_vec(M.x)];
 
% add [0] states if not specified
%--------------------------------------------------------------------------
try
    M.f = spm_funcheck(M.f);
catch
    M.f = @(x,u,P,M) sparse(0,1);
    M.x = sparse(0,0);
end

 
% output nonlinearity, if specified
%--------------------------------------------------------------------------
try
    g   = spm_funcheck(M.g);
catch
    g   = @(x,u,P,M) x;
    M.g = g;
end
 
% % Bilinear approximation (1st order)
% %--------------------------------------------------------------------------
% [M0,M1] = spm_bireduce(M,P);
% m       = length(M1);                     % m inputs
 
% delays
%--------------------------------------------------------------------------
try
    D  = max(round(M.delays/U.dt),1);
catch
    D  = ones(M.l,1)*round(u/v);
end


% Evaluation times (t) and indicator array for inputs (su) and output (sy)
%==========================================================================
 
% get times that the input changes
%--------------------------------------------------------------------------
%i     = [1 (1 + find(any(diff(U.u),2))')];
su    = ones(1,u);%sparse(1,i,1,1,u);
 
% get times that the response is sampled
%--------------------------------------------------------------------------
s     = ceil((0:v - 1)*u/v);
for j = 1:M.l
    i       = s + D(j);
    sy(j,:) = sparse(1,i,1:v,1,u);
end
 
% time in seconds
%--------------------------------------------------------------------------
t     = find(su | any(sy));
su    = full(su(:,t));
sy    = full(sy(:,t));
dt    = [diff(t) 0]*U.dt;
 
 
% Integrate
%--------------------------------------------------------------------------
y     = zeros(M.l,v);
X     = zeros(length(x),v);
%J     = M0;
U.u   = full(U.u);

for i = 1:length(t)

      q          = spm_unvec(x,M.x);
     [fx dfdx]   = M.f(q,U.u(t(i),:),P,M);
    % dfdx        = spm_cat(spm_diff(M.f,x,U.u(t(i),:),P,M,1));

    % output sampled
    %----------------------------------------------------------------------
    if any(sy(:,i))
        
        j      = find(sy(:,i));
        s      = sy(j(1),i);
        U.asl   = M.asl*(-1)^(s);
        q      = spm_vec(M.g(q,U.asl,P,M));
        y(j,s) = q(j);
    end
 
    % compute updated states x = expm(J*dt)*x;
    %----------------------------------------------------------------------
    %x  = spm_expm(J*dt(i),x);
    
        % update dx = (expm(dt*J) - I)*inv(J)*fx
    %----------------------------------------------------------------------
   % x      = spm_unvec(spm_vec(x) + spm_dx(dfdx,fx(:),dt(i)),x);
   % Ito-Taylor
    x         = x + dt(i)*fx(:) + 0.5*dt(i)^2*(full(dfdx)*fx(:));

    
    % check for convergence
    %----------------------------------------------------------------------
    if norm(x,1) > 1e6, break, end
    X(:,i) = x(:);
    if nargout>2
      dFx(i,:) = fx;   
    end
 
end
y      = real(y');
X      = real(X');
% X(:,3:5) = exp(X(:,3:5));
% figure(100), subplot(211),plot(X(:,3:5)); drawnow;
%              subplot(212),plot(X(:,[1 2 6])); drawnow; 
