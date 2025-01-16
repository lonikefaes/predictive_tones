function [M] = LBR_param_priors(N,K,A,B,C),
% LBR_parameters defines structure of parameters for laminar BOLD response 
%                (LBR) model (see LBR_model.m). It applies parameter values
%                 as describe in Table 1 of Havlicek, M. & Uludag, K. (2019) BioRxiv  
%
% INPUT:  K - Number of cortical depths
%
% OUTPUT: P - structure with all default parameters for LBR model
%
% AUTHOR: Martin Havlicek, 5 August, 2019
%
% REFERENCE: Havlicek, M. & Uludag, K. (2019) A dynamical model of the
%            laminar BOLD response, BioRxiv, doi: https://doi.org/10.1101/609099 
%
% EXAMPLE:
%            K = 6;
%            P = LBR_parameters(K);
%            disp(P);
%--------------------------------------------------------------------------
M.N  = N;    % Number of depths
M.K  = K;

% Neuronal parameter:
%--------------------------------------------------------------------------
P0.sigma = 3;
P0.mu    = 1.5;
P0.lam   = 0.1;
P0.A     = A;
P0.B     = B;
P0.C     = C;
P0.Bmu   = 0;
P0.Blam  = 0;

% NVC parameters:
% --------------------------------------------------------------------------
P0.c1      = 0.6;
P0.c2      = 1.5;
P0.c3      = 0.6;  

depths     = linspace(0,100,2*K+1); % Normalized distance to the center of individual depths (in %)
M.depths   = depths(2:2:end)';
P0.nsig    = 0.005;

% LAMINAR HEMODYNAMIC MODEL:
%--------------------------------------------------------------------------
% Baseline physiological parameters:
P0.V0t   = 3.5;  % Total (regional) amount of CBV0 in the gray matter (in mL) [1-6]
P0.w_v = 0.5;  % CBV0 fraction of microvasculature (i.e. venules here )with respect to the total amount 
P0.x_v = [];   % CBV0 fraction across depths in venules 
P0.x_d = [];   % CBV0 fraction across depths in ascending veins
P0.s_v = 0;    % Slope of CBV increase (decrease) in venules [0-0.3]
P0.s_d = 0.4;  % Slope of CBV increase in ascending vein     [0-1.5]
P0.s_d0 = 0;  % Slope of CBV increase in ascending vein     [0-1.5]
P0.s_d2 = 0;  % Slope of CBV increase in ascending vein     [0-1.5]

P0.t0v = 1;    % Transit time through microvasculature(in second)
P0.E0v = 0.35; % Baseline oxygen extraction fraction in venules
P0.E0d = 0.35; % Baseline oxygen extraction fraction in venules
P0.E0p = 0.35; % Baseline oxygen extraction fraction in venules

% Parameters describing relative relationship between physiological variable:
% CBF-CBV coupling (steady-state)

P0.al_v = 0.3; % For venules
P0.al_d = 0.2; % For ascending vein


% CBF-CMRO2 coupling (steady-state)
P0.nr = 4;         % n-ratio   (Ref. Buxton et al. (2004) NeuroImage)

% CBF-CBV dynamic uncoupling 
P0.tau_v_in = 2; % For venules - inflation 
P0.tau_v_de = 2; %             - deflation
P0.tau_d_in = 2; % For AV      - inflation 
P0.tau_d_de = 2; %             - deflation

% LAMINAR BOLD SIGNAL MODEL:
%--------------------------------------------------------------------------
M.TE     = 0.028;     % echo-time (in sec)

% Hematocrit fraction
P0.Hct_v  = 0.35;      % For venules, Ref. Lu et al. (2002) NeuroImage
P0.Hct_d  = 0.38;      % For ascending vein
P0.Hct_p  = 0.42;      % For pial vein


M.B0      = 7;                 % Magnetic field strenght (in Tesla)  
P0.gyro   = 2*pi*42.6*10^6;    % Gyromagnetic constant for Hydrogen
P0.suscep = 0.264*10^-6;       % Susceptibility difference between fully oxygenated and deoxygenated blood

% Water proton density:
P0.rho_t  = 0.89;                % For gray matter tissue 
P0.rho_v  = 0.95 - P0.Hct_v*0.22; % For blood (venules) Ref. Lu et al. (2002) NeuroImage
P0.rho_d  = 0.95 - P0.Hct_d*0.22; % For blood (ascending vein)
P0.rho_p  = 0.95 - P0.Hct_p*0.22; % For blood (pial vein)
P0.rho_tp = 0.95;                % For gray matter tissue % CSF   

% Relaxation rates for 7 T (in sec-1)
P0.R2s_t  = 34; % For gray matter tissue
P0.R2s_v  = 80; % For blood (venules)
P0.R2s_d  = 85; % For blood (ascending vein)
P0.R2s_p  = 90; % For blood (pial vein)

% Slope of change in R2* of blood with change in extraction fration during activation 
P0.r0v    = 228;      % For venules   
P0.r0d    = 232;      % For ascending vein
P0.M0     = 100;

M.P0 = P0;
M.x  = zeros(N*4+K*4,1);
M.xn = zeros(N,4);
M.xk = zeros(K,4);