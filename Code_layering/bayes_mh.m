function [w C logEV] = bayes_reg(yd,X,mode),

p = size(X,2);
N = size(X,1);
% measurement noise covariance (gamma prior)
beta0  = 1e-2;
alpha0 = 1e-2;
alpha  = alpha0;
beta   = beta0;

R      = eye(N)*(beta0/alpha0); 
% Shrinkage priors
wp     = zeros(p,1);         % zero mean
Cp     = eye(p)*1;           % covariance  


% posterior:
w     = zeros(p,1);   % mean
C     = zeros(p,p);   % covariance


if mode
    ITER = 20;
else
    ITER = 20;
end

% Loop
for it = 1:ITER
    
    % E-step:
    %----------------------------------------------------------------------
    C = inv(X'*inv(R)*X + inv(Cp));
    w = C*(X'*inv(R)*yd + inv(Cp)*wp);
    
    % M-step(s):
    %----------------------------------------------------------------------
    if ~mode
        % update for noise covariance:
        gamma = N - alpha*trace(C);
        alpha = p/((w'*w) + trace(C));
        beta  = (N - gamma)/((yd - X*w)'*(yd - X*w));
        R     = eye(N)/beta;
        % update for prior covariance (one parameter)
        Cp    = eye(N)/alpha;
    else % OR...
        % update for noise covariance:
        alpha = alpha + N/2;
        beta  = beta  + (yd - X*w)'*(yd - X*w)/2;
        R     = eye(N)*(beta/alpha);
        % update for prior covariance (simple ARD):
        Cp    = diag(diag(w*w' + C));
    end
    
    
%     figure(f1);
%     subplot(311); plot(X*w,'g');
%     axis([0 101 -1 2])
%     subplot(312); plot(X*diag(w));
%     b = X*diag(w);
%     axis([0 101 min(b(:))-0.01 max(b(:))+0.01])
% 
%     subplot(313); 
%     if ~mode,  
%         bar(diag(Cp));  
%         axis([0 101 0 max(diag(Cp)+0.01)])
%     else, 
%         bar(diag(Cp)/sum(diag(Cp))); 
%         axis([0 101 0 1])
%     end;
   % drawnow;

    e1  = yd  - X*w;
    e2  = wp  - w;
    p   = length(w);
    
    logdetR	   = 2*sum(log(diag(R)));   % approximation to the determinant (stable)...
    logdetC	   = 2*sum(log(diag(C)));
    logdetCp   = 2*sum(log(diag(Cp)));
    accuracy   = -0.5*N*logdetR-0.5*(e1'/R*e1);
    complexity =  0.5*p*logdetCp-0.5*logdetC+0.5*(e2'/(Cp)*e2);
    const      = -0.5*N*log(2*pi);
    
    logEV      =  accuracy - complexity + const;
    logLik     =  accuracy + const;
   % disp(logEV(it));
%    figure(f2), plot(logEV); ylabel('log-Evidence'); xlabel('Iteration');

  
end