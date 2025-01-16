function [f dfdx] = LBR_model_fx(x,u,P,M)
% LBR_model calculates laminar BOLD response (LBR) as in detail described by 
%           Havlicek, M. & Uludag, K. (2019) BioRxiv  
%
% INPUTS:
%       P - structure of model parameters (see LBR_parameters.m)
%
%       cbf - matrix defining laminar cerebral blood flow (CBF) response, (time,depths). (Required);
%
%       cmro2 - matrix defining laminar changes in oxygen metabolism (CMRO2), (time,depth). 
% OUTPUTS:
%       LBR - matrix containing laminar BOLD responses in perecent signal change (time,depths)
%   
%       Y - structure with all baseline and relative physiological
%       variables underlying BOLD response
%
%       LBRpial - BOLD response of the pial vein in perecent signal change (0th depth) (time,1)
%
% AUTHOR: Martin Havlicek, 5 August, 2019
%
% REFERENCE: Havlicek, M. & Uludag, K. (2019) A dynamical model of the
%            laminar BOLD response, BioRxiv, doi: https://doi.org/10.1101/609099 
%
% EXAMPLE:
%        For steady-state:
%               K = 6;                       % Number of depths
%               P = LBR_parameters(K);       % Get parameter structure with default values
%               cbf = ones(P.T/P.dt,K)*1.6;  % Define model input (Relative blood flow across depths)
%               P.s_d = 0.4;                 % Define the slope of increase of CBV0 in the ascening vein
%               [LBR,Y] = LBR_model(P,cbf);  % Generate the LBR
%               figure(1), 
%               plot(P.l,flipud(LBR(end,:)')); % plot the LBR profile as a function of normalized cortical depth
%               xlim([0 100]); ylim([0 6]); xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
%
%        For dynamic response:
%               K = 6;  
%               P.N = neuronal_NVC_parameters(K);  % consider default parameters
%               P.N.T = 30; % (in seconds)
%               P.H = LBR_parameters(K);
%               P.H.T  = P.N.T;
%               P.H.dt = P.N.dt;
%               U.u = zeros(P.N.T/P.N.dt,K);
%               dur = 2/P.N.dt; % 2 sec stimulus
%               onset     = 2/P.N.dt;
%               offset    = onset + dur;
%               U.u(onset:offset,:) = 1;
%               [neuro, cbf]  = neuronal_NVC_model(P.N,U);
%               P.H.alpha_v   = 0.35;
%               P.H.alpha_d   = 0.2;
%               P.H.tau_d_de  = 30;
%               [LBR,Y]       = LBR_model(P.H,cbf);
%               time_axis = [0:P.H.dt:P.H.T-P.H.dt];
%               figure(1), 
%               subplot(121), plot(time_axis,cbf); xlim([time_axis(1), time_axis(end)]); ylim([0.5 2]); 
%                            xlabel('1 - Cortical depth (%)'); ylabel('Relative CBF (%)'); axis square;
%               subplot(122), plot(time_axis,LBR); xlim([time_axis(1), time_axis(end)]); ylim([-1 4]);  %                           xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
%
%--------------------------------------------------------------------------

persistent fx,

%%
% Neuronal model parameters:
%--------------------------------------------------------------------------

N       = M.N;      % Number of neuronal depths 
P0      = M.P0;

A       = P.A;
C       = P.C;
B       = P.B;
Bmu     = P.Bmu;
Blam    = P.Blam;

sigma = P0.sigma*exp(P.sigma);
A     = A - diag(diag(A)) - diag(sigma*exp(diag(A)));

for i = 1:size(B,3)
    A = A + u(i)*B(:,:,i);
end

mu    = P.mu;
lam   = P.lam;
for i = 1:size(Bmu,2), 
    mu = mu + Bmu(:,i)*u(i);
end
if length(mu) == N,
    mu = P0.mu.*exp(mu);
else
    mu = P0.mu*ones(N,1)*exp(mu);
end

for i = 1:size(Blam,2),         % different rows are associated with different inputs
    lam = lam + Blam(:,i)*u(i); % same inputs for mu and lambda
end
if length(lam) == N,
    lam = P0.lam.*exp(lam);
else
    lam = P0.lam*ones(N,1)*exp(lam);
end
CU  = C*u(:);


% NVC parameters:
% --------------------------------------------------------------------------
if length(P.c1) == N,
    c1 = P0.c1.*exp(P.c1);
else
    c1 = P0.c1*ones(N,1)*exp(P.c1); 
end
if length(P.c2) == N,
    c2 = P0.c2.*exp(P.c2);
else
    c2 = P0.c2*ones(N,1)*exp(P.c2); 
end
if length(P.c3) == N,
    c3 = P0.c3.*exp(P.c3);
else
    c3 = P0.c3*ones(N,1)*exp(P.c3); 
end

%%
% Hemodynamic model parameters:
%--------------------------------------------------------------------------
K      = M.K;      % Number of depths (Reference: Superfical depth close to CSF is k = 1;)
depths = M.depths;

if N==K
    n2k = eye(N);
else
    m1      = linspace(0,1,2*N+1);
    m1      = m1(2:2:end);
    m2      = linspace(0,1,N+2);
    m2      = m2(2:end-1);
    m       = (m1+m2)/2;
    dm     = diff(m);
    
    m      = repmat([m(1)-dm(1),m,m(end)+dm(1)],K+2,1);
    if P.s~=0
        m   = m + 0.01*(m-0.5)*exp(P.s);
    end
    
    nsig   = ((1/(3*(N+1)))^2);   % sigma = 
    nsig   = ((nsig*1.2 - nsig*0.8)*normcdf(P.nsig) + nsig*0.8).*ones(K+2,N+2);
    sp     = linspace(0,1,2*K+1);
    sp     = sp(2:2:end);
     dsp     = diff(sp);
    sp      = repmat([sp(1)-dsp(1),sp,sp(end)+dsp(1)]',1,N+2);

       
    

    n2k    = 1./sqrt(2*pi*nsig).*exp(-(sp-m).^2./(2.*nsig));
    n2k    = n2k./repmat(sum(n2k,2),1,size(n2k,2));
   % n2k    = [max(n2k(:,1),n2k(:,2)),n2k(:,3:end-2),max(n2k(:,end-1),n2k(:,end))];
    n2k    = n2k(2:end-1,2:end-1);
    n2k(1,1)= n2k(1,1) + 0.05*exp(P.nb);
   % n2k1    = n2k
    n2k(end,end)= n2k(end,end) + 0.05*exp(P.nb);
   % n2k    = n2k./repmat(sum(n2k,2),1,N);
  %  n2k    = n2k./repmat(max(n2k),K,1);
end

% BASELINE PARAMETERS:
V0t       = P0.V0t.*exp(P.V0t);    % Total amount of CBV0 within GM tisse (in mL)
%V0t_p     = P.V0t_p;  % Total amount of CBV0 in pial vein (in mL)

w_v       = P0.w_v.*exp(P.w_v);    % Fraction of CBV0 in venules with respect to the total, 
w_d       = 1-w_v;    % Fraction of CBV0 in ascending vein with respect to the total, 

s_v       = P0.s_v.*exp(P.s_v);    % Slope of CBV0 increase towards the surface in venules 
%s_d       = (0.8 - 0.0)*normcdf(P.s_d);

s_d       = P0.s_d.*exp(P.s_d);    % Slope of CBV0 increase towards the surface in ascending vein 
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
    x_d  = 10+s_d*flipud(depths(:)) + s_d2*max(0,flipud(depths(:))-s_d0);  % Possibility to define linear increase 
end  
x_d      = x_d./sum(x_d);          % Fraction of CBV0 across depths in venules 

V0v      = V0t*w_v*x_v;            % CBV0 in venules
V0d      = V0t*w_d*x_d;            % CBV0 in ascending vein
%V0p      = V0t_p;                  % CBV0 in pial vein

% Transit time through venules (or microvasculature in general)
if length(P.t0v) == K,
    t0v    = P0.t0v.*exp(P.t0v);              % depth-specific defined by user
else
    t0v    = P0.t0v*ones(K,1).*exp(P.t0v);   % default
end;

% Depth-specific baseline CBF:  
F0v     = V0v./t0v;          % Note: can be also defined directly and t0v calculated from V0v and F0v
F0d     = flipud(cumsum(flipud(F0v)));

% Depth-specific transit time:  
t0v     = V0v./F0v;   
t0d     = V0d./F0d;


%% PARAMETERS DESCRIBING RELATIVE RELATIONSHIPS BETWEEN PHYSIOLOGICAL VARIABLES:
%
% n-ratio (= (cbf-1)./(cmro2-1)). Not used if cmro2 response is directly specified as an input
if length(P.nr) == K,             % For venules (microvasculature)
    nr      = P0.nr.*exp(P.nr);                % Depth-specific defined by user
else
    nr      = P0.nr*ones(K,1)*exp(P.nr);      % Default
end;

% Grubb's exponent alpha (i.e CBF-CBV steady-state relationship)
if length(P.al_v) == K,       % For venules
    al_v    = P0.al_v.*exp(P.al_v);             % Depth-specific defined by user 
else
    al_v    = P0.al_v*ones(K,1).*exp(P.al_v);  % Default
end;
if length(P.al_d) == K,       % For ascending vein
    al_d    = P0.al_d.*exp(P.al_d);             % Depth-specific defined by user  
else
    al_d    = P0.al_d*ones(K,1).*exp(P.al_d);  % Default
end

% CBF-CBV uncoupling (tau) during inflation and deflation:
if length(P.tau_v_in) == K,      % For venules (inflation)
    tau_v_in  = P0.tau_v_in.*exp(P.tau_v_in);             % Depth-specific defined by user
else
    tau_v_in  = P0.tau_v_in*ones(K,1).*exp(P.tau_v_in);  % Default
end
if P0.tau_v_same == 0
    if length(P.tau_v_de) == K,      % For venules (inflation)
        tau_v_de  = P0.tau_v_de.*exp(P.tau_v_de);             % Depth-specific defined by user
    else
        tau_v_de  = P0.tau_v_de*ones(K,1).*exp(P.tau_v_de);  % Default
    end
else
    if length(P.tau_v_de) == K,      % For venules (inflation)
        tau_v_de  = P0.tau_v_in.*exp(P.tau_v_in);             % Depth-specific defined by user
    else
        tau_v_de  = P0.tau_v_in*ones(K,1).*exp(P.tau_v_in);  % Default
    end
end

if length(P.tau_d_in) == K,      % For ascending vein (inflation)
    tau_d_in  = P0.tau_d_in.*exp(P.tau_d_in);             % Depth-specific defined by user 
else
    tau_d_in  = P0.tau_d_in*ones(K,1)*exp(P.tau_d_in);   % Default  
end;
if P0.tau_d_same == 0
    if length(P.tau_d_de) == K,      % For ascending vein (inflation)
        tau_d_de  = P0.tau_d_de.*exp(P.tau_d_de);             % Depth-specific defined by user
    else
        tau_d_de  = P0.tau_d_de*ones(K,1)*exp(P.tau_d_de);   % Default
    end;
else
    if length(P.tau_d_de) == K,      % For ascending vein (inflation)
        tau_d_de  = P0.tau_d_in.*exp(P.tau_d_in);             % Depth-specific defined by user
    else
        tau_d_de  = P0.tau_d_in*ones(K,1)*exp(P.tau_d_in);   % Default
    end;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tau_v = tau_v_in;
tau_d = tau_d_in;


if ~isempty(fx)
    dv_v  = fx(4*N+1:4*N+K);
    dv_d  = fx(4*N+2*K+1:4*N+3*K);
    tau_v(dv_v<0) = tau_v_de(dv_v<0);
    tau_d(dv_d<0) = tau_d_de(dv_d<0);
end


xn      = spm_unvec(x(1:N*4),M.xn);
xk      = spm_unvec(x(N*4+1:end),M.xk);
xn(:,4) = exp(xn(:,4));
xk      = exp(xk);



IN      = num2cell([xn(:);...
                    xk(:);...
                    A(:);...
                    mu(:);...
                    CU(:);...
                    lam(:);...
                    c1(:);...
                    c2(:);...
                    c3(:);...
                    n2k(:);...
                    al_v(:);...
                    al_d(:);
                    t0v(:);...
                    t0d(:);...
                    nr(:);...
                    tau_v(:);...
                    tau_d(:);...
                    F0v(:);...
                    F0d(:)]);
                  if N==3,
                    if K>=3 || K<=15
                        f       = eval(['sym_fx_LBR_N3K',num2str(K),'(IN{:})']);
                        dfdx    = eval(['sym_dfdx_LBR_N3K',num2str(K),'(IN{:})']);
                    else
                        disp('model not defined');
                    end
                elseif N == 6
                       if K>=6 || K<=15
                        f       = eval(['sym_fx_LBR_N6K',num2str(K),'(IN{:})']);
                        dfdx    = eval(['sym_dfdx_LBR_N6K',num2str(K),'(IN{:})']);   
                       else
                        disp('model not defined');  
                      end
                elseif N == 2
                       if K>=6 || K<=12
                        f       = eval(['sym_fx_LBR_N2K',num2str(K),'(IN{:})']);
                        dfdx    = eval(['sym_dfdx_LBR_N2K',num2str(K),'(IN{:})']);   
                       else
                        disp('model not defined');  
                      end
       
                else
                    disp('model not defined')
                end

                fx = f;

