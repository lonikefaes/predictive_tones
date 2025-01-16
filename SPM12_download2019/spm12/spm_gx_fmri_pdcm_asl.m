function [g,dgdx] = spm_gx_fmri_pdcm_asl(x,u,P,M)
% Simulated BOLD response to input
% FORMAT [g,dgdx] = spm_gx_fmri(x,u,P,M)
% g          - BOLD response (%)
% x          - state vector     (see spm_fx_fmri)
% P          - Parameter vector (see spm_fx_fmri)
% M          - model specification structure (see spm_nlsi)
%__________________________________________________________________________
%
% This function implements the BOLD signal model described in: 
%
% Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
% Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston & Klaas Enno Stephan
% $Id: spm_gx_fmri.m 6262 2014-11-17 13:47:56Z karl $
 
 
% Biophysical constants for 1.5T
%==========================================================================
 
% time to echo (TE) (default 0.04 sec)
%--------------------------------------------------------------------------
n = size(x,1);

if M.asl~=0 && u~=0
    TE  = M.TEasl;%.*exp(P.TEasl);
else
    TE  = M.TE;
end
% resting venous volume (%)
%--------------------------------------------------------------------------
if M.asl~=0 && u~=0
    if length(P.V0asl)<n
        V0  = 0.04*ones(n,1)*exp(P.V0asl);
    else
        V0  = 0.04*exp(P.V0asl);
    end
else
     if length(P.V0)<n
        V0  = 0.04*ones(n,1)*exp(P.V0);
    else
        V0  = 0.04*exp(P.V0);
    end   
    
end
% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation S:  R_iv = r0*[(1 - S)-(1 - S0)] (Hz)
%--------------------------------------------------------------------------

 
% resting oxygen extraction fraction
%--------------------------------------------------------------------------
E0  = 0.35;
     % Susceptibility difference between fully oxygenated and deoxygenated blood


Hct  = 0.38;       % Hematocrit fraction

B0     = M.B0;              % Field strenght        
gyro   = 2*pi*42.6*10^6;    % Gyromagnetic constant 
suscep = 0.264*10^-6;       % Susceptibility difference

nu0   = suscep*gyro*Hct*B0;


% Water proton density 
rho_t  = 0.89;  % In GM tissue
rho_b  = 0.95 - Hct*0.22;  % In blood  Ref. Lu et al. (2002) NeuroImage


% Relaxation rates (in sec-1):

R2s_t  = 34;         % For tissue (for 7T)


R2s_b  = 220;           % For venous blood (for 7T)

% (Baseline) Intra-to-extra-vascular signal ratio
ep   = rho_b./rho_t.*exp(-TE*R2s_b)./exp(-TE*R2s_t);        % For venules

% Slope of change in R2* of blood with change in extraction fration during activation 
r0    = 228;      % For 7T

%-Coefficients in BOLD signal model
%==========================================================================

k1     = 4.3.*nu0.*E0.*TE;
k2     = ep.*r0.*E0.*TE;
k3     = 1 - ep; 
%k3     = ep - 1; 

%-Output equation of BOLD signal model
%==========================================================================
f   = exp(x(:,3));
v   = exp(x(:,4));
q   = exp(x(:,5));
bold  = V0.*(k1.*(1 - q) + k2.*(1 - q./v) + k3.*(1-v));
if M.asl~=0 && u~=0
%     T1a    = 1.664;
%     T1t    = 1.122;
%     TI     = 2.35;
   % alpha = 0.98;
    CBF0 =  M.ASLb.*exp(P.fbase);
    M0   = M.M0asl.*exp(P.M0asl);
  %  g  = (bold+1).*M0 + u.*f.*CBF0/2;
%      lambda = 0.9;
%      T1app  = 1./T1t + f.*CBF0/lambda;
     % g = M0.*(1-exp(-TR./T1app)).*(1+bold);
      if u == 1
%          %fmod = ;  % tag control labeling (-1 +1)
% %         %
%         
          g = M0.*(1+bold);
      elseif u == -1;
%          g = M0.*(1-exp(-TR./T1app)...
%                - 2*f.*CBF0/lambda.*(exp(-TR./T1app)-exp(-TI./T1a))./(1./T1a - 1./T1app)).*(1+bold);
%             % - 2*f.*CBF0/lambda.*T1app.*(1-exp(-TR./T1app)).*exp(-TI./T1a)).*(1+bold);
          g = M0.*(1 - 2*f.*CBF0).*(1+bold);  
      end
    
%     SnB = CBF0.*f.*(1 - alpha*(1 + u)).*exp(-TI./T1a);
%     SnS = M.M0asl.*exp(P.M0asl);
%     g   = (SnB+SnS).*(bold+1);
    
else
    g  = (bold+1).*M.M0.*exp(P.M0);
end
%asl  = boldc;% + fmod;

% if nargout == 1, return, end
% 
% 
% %-derivative dgdx
% %==========================================================================
% [n m]      = size(x);
% dgdx       = cell(1,m);
% [dgdx{:}]  = deal(sparse(n,n));
% dgdx{1,4}  = diag(-V0*(k3.*v - k2.*q./v));
% dgdx{1,5}  = diag(-V0*(k1.*q + k2.*q./v));
% dgdx       = spm_cat(dgdx);
