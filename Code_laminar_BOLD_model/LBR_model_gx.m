function [LBR] = LBR_model_gx(x,u,P,M)

K  = M.K;
N  = M.N;
P0 = M.P0;
depths = M.depths;
%% Parameters for laminar BOLD signal equation (for 7 T field strenght):
%--------------------------------------------------------------------------
% Baseline CBV in fraction with respect to GM tissue

V0t       = P0.V0t.*exp(P.V0t);    % Total amount of CBV0 within GM tisse (in mL)

w_v       = P0.w_v.*exp(P.w_v);    % Fraction of CBV0 in venules with respect to the total,
w_d       = 1-w_v;    % Fraction of CBV0 in ascending vein with respect to the total,

s_v       = P0.s_v.*exp(P.s_v);    % Slope of CBV0 increase towards the surface in venules
s_d       = P0.s_d.*exp(P.s_d);    % Slope of CBV0 increase towards the surface in ascending vein
%s_d       = (0.8 - 0.0)*normcdf(P.s_d);
s_d2      = P0.s_d2.*exp(P.s_d2);
s_d0      = P0.s_d0.*exp(P.s_d0);
% Depth-specific CBV0:
if length(P.x_v) == K,             % For venules
    x_v  = P0.x_d.*exp(P.x_v);                  % Depth-specific fractions defined by user
else
    x_v  = 10+s_v*flipud(depths(:));  % Possibility to define linear increase (default s_v = 0)
end
x_v      = x_v./sum(x_v);          % Fraction of CBV0 across depths in venules

if length(P.x_v) == K,             % For ascending vein
    x_d  = P0.x_d.*exp(P.x_d);                  % Depth-specific fractions defined by user
else
    x_d  = 10+s_d*flipud(depths(:)) + s_d2*(flipud(depths(:)-s_d0));  % Possibility to define linear increase
end
x_d      = x_d./sum(x_d);          % Fraction of CBV0 across depths in venules

V0v      = V0t*w_v*x_v;            % CBV0 in venules
V0d      = V0t*w_d*x_d;            % CBV0 in ascending vein


V0vq = V0v./100*K;
V0dq = V0d./100*K;


% Baseline oxygen extraction fraction
if length(P.E0v) == K,
    E0v        = P0.E0v.*exp(P.E0v);     % depth-specific defined by user
else
    E0v        = P0.E0v*ones(K,1).*exp(P.E0v);     % default
end
if length(P.E0d) == K,
    E0d        = P0.E0d.*exp(P.E0d);     % depth-specific defined by user
else
    E0d        = P0.E0d*ones(K,1).*exp(P.E0d);    % default
end


TE     = M.TE;          % echo-time (sec)

Hct_v  = P0.Hct_v;       % Hematocrit fraction
Hct_d  = P0.Hct_d;

B0     = M.B0;          % Field strenght
gyro   = P0.gyro;        % Gyromagnetic constant
suscep = P0.suscep;      % Susceptibility difference

nu0v   = suscep*gyro*Hct_v*B0;
nu0d   = suscep*gyro*Hct_d*B0;


% Water proton density
rho_t  = P0.rho_t;  % In GM tissue
rho_v  = P0.rho_v;  % In blood (venules) Ref. Lu et al. (2002) NeuroImage
rho_d  = P0.rho_d;  % In blood (ascening vein)

% Relaxation rates (in sec-1):
if length(P0.R2s_t) == K,  % For tissue
    R2s_t  = P0.R2s_t;           %
else
    R2s_t  = ones(K,1).*P0.R2s_t;   % (sec-1)
end
if length(P0.R2s_v) == K,  % For venules
    R2s_v  = P0.R2s_v;               % (sec-1)
else
    R2s_v  = ones(K,1)*P0.R2s_v; % (sec-1)
end
if length(P0.R2s_d) == K,  % For ascening vein
    R2s_d  = P0.R2s_d;           % (sec-1)
else
    R2s_d  = ones(K,1)*P0.R2s_d; % (sec-1)
end

% (Baseline) Intra-to-extra-vascular signal ratio
ep_v   = rho_v./rho_t.*exp(-TE*R2s_v)./exp(-TE*R2s_t);        % For venules
ep_d   = rho_d./rho_t.*exp(-TE*R2s_d)./exp(-TE*R2s_t);        % For ascending vein


% Slope of change in R2* of blood with change in extraction fration during activation
r0v    = P0.r0v;      % For venules
r0d    = P0.r0d;      % For ascending vein


H0     = 1./(1 - V0vq - V0dq + ep_v.*V0vq + ep_d.*V0dq);  % constant in front

% H0v    = 1./(1 - V0vq + ep_v.*V0vq);
% H0d    = 1./(1 - V0dq + ep_d.*V0dq);


k1v     = 4.3.*nu0v.*E0v.*TE;
k2v     = ep_v.*r0v.*E0v.*TE;
k3v     = 1 - ep_v;

k1d     = 4.3.*nu0d.*E0d.*TE;
k2d     = ep_v.*r0d.*E0d.*TE;
k3d     = 1 - ep_d;

M.xk(:) = x(N*4+1:end);
M.xk    = exp(M.xk);

v_v     = M.xk(:,1);
v_d     = M.xk(:,3);
q_v     = M.xk(:,2);
q_d     = M.xk(:,4);

LBR = H0.*((1-V0vq-V0dq).*(k1v.*V0vq.*(1-q_v)       +k1d.*V0dq.*(1-q_d))+...
    +k2v.*V0vq.*(1-q_v./v_v) +k2d.*V0dq.*(1-q_d./v_d)+...
    +k3v.*V0vq.*(1-v_v)      +k3d.*V0dq.*(1-v_d)).*100;

if N~=K
    kernel = M.kernel;
    %LBRi = interp1([1:K]',flipud(LBR),[0:K+1]','linear','extrap');
    LBR2 = conv([LBR(end);flipud(LBR);LBR(1)],kernel(:),'same');
    LBR2 = flipud(LBR2(2:end-1));
    LBR3 = conv([LBR(1);LBR;LBR(end)],kernel(:),'same');
    LBR3 = LBR3(2:end-1);
    LBR  = max(LBR2,LBR3);
end

%LBR2 = conv(flipud(LBR),kernel(:),'same');       
%LBR = LBR.*P0.M0.*exp(P.M0);  
% LBRv = H0v.*((1-V0vq).*(k1v.*V0vq.*(1-q_v)+k2v.*V0vq.*(1-q_v./v_v)+...
%                            +k3v.*V0vq.*(1-v_v))).*100;
%                        
% LBRd = H0d.*((1-V0dq).*(k1d.*V0dq.*(1-q_d))+k2d.*V0dq.*(1-q_d./v_d)+...
%                            +k3d.*V0dq.*(1-v_d)).*100;

