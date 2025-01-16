function [BOLD, CBF,NEU, BOLD_venules, BOLD_drain, Vv, Vd,dt] = laminar_bold_model(M,Bin,Bmu,Bsig)
  
%clear all; close all;

K         = M.K;      % Number of depths
% Reference: Superfical depth close to CSF is K = 1;

V0t       = 3.5;      % Total baseline CBV (withi GM)
V0_p      = 4;        % baseline CBV pail
w_v       = 0.5;      % Volume fraction microvasculatur [venules+capillaries], 
w_d       = 1-w_v;     % Volume fraction draining vein;
s_d       = 0.5;      % V0 slope increase in draining vein 
r_v       = 1;        % V0 slope increase in microvasculature (defined as ratio superficial/deep)   

% Note: Can be defined also as depth specific (but assumed constant for now)
n          = ones(K,1)*4;    % n-ratio (rCBF-1)/(rCMRO2-1)  
tt_v       = ones(K,1)*1;    % transit time through microvasculature
E0v         = ones(K,1)*0.35; % baseline oxygen extraction fraction
E0d         = ones(K,1)*0.35; % baseline oxygen extraction fraction
E0p         = ones(1,1)*0.35; % baseline oxygen extraction fraction


% Grubb's alpha (CBF-CBV steady-state relationship)
alpha_v    = ones(K,1)*0.35;
alpha_d    = ones(K,1)*0.25;
alpha_p    = ones(1,1)*0;

% CBF-CBV uncoupling (tau) during inflation and deflation:
tau_v_in  = ones(K,1)*10;  % microvasculature
tau_v_de  = ones(K,1)*10; 
tau_d_in  = ones(K,1)*20;  % draining vein
tau_d_de  = ones(K,1)*40; 
tau_p_in  = ones(1,1)*40;  % pial vein
tau_p_de  = ones(1,1)*100;

Simf     = [1.6];

% Hemodynamic model parameters:
%--------------------------------------------------------------------------
% baseline CBV:
V0v      = V0t*w_v;%(ml/100g)
V0d      = V0t*w_d;
V0p      = V0_p;

trend_v  = linspace(r_v,1,K)';
V0vq     = V0v.*((trend_v./sum(trend_v)))/100*K;

trend_d  = flipud(cumsum(flipud([(s_d).*V0vq(1:K-1);V0vq(K)])));
V0dq     = V0d.*(trend_d./sum(trend_d))./100*K;
V0pq     = V0p./100;   

% baseline CBF:  
F0vq     = V0vq./tt_v; % Note: can be also defined directly 
F0dq     = flipud(cumsum(flipud(F0vq)));
F0pq     = F0dq(1);

% transit times (below expressed directly as ratios):
tt_v     = V0vq./F0vq;   
tt_d     = V0dq./F0dq;
tt_q     = V0pq./F0pq;


tau_v = tau_v_in;
tau_d = tau_d_in;
tau_p = tau_p_in;

% BOLD signal equation parameters:
%--------------------------------------------------------------------------
TE     = 0.028;     % echo-time (sec)
r0v    = 128;   % ???recalculate     % ~for 7 T
r0d    = 132;  % A* + B* (1-Y) + C* (1-Y)^2; Hct 21: A*=14.5; B*=33.2; C*=55.4 (at 3T)
r0p    = 136;  % Hct 44: A*=17.5; B*=39.1; C*=119
Hct_v  = 0.35; % Hct 57: A*=18.3; B*=56.0; C*=128
Hct_d  = 0.38;
Hct_p  = 0.42;
B0     = 7;
gyro   = 2*pi*42.6*10^6;
suscep = 0.264*10^-6;
nu0v   = suscep*gyro*Hct_v*B0;
nu0d   = suscep*gyro*Hct_d*B0;
nu0p   = suscep*gyro*Hct_p*B0; % for 7 T

ep_v   = 0.1;       % intra-to-extravascular baseline signal ratio for draining veins
ep_d   = 0.2;       % intra-to-extravascular baseline signal ratio for draining veins
ep_p   = 0.3;
Vi0    = V0vq + V0dq;  % total intravascular baseline blood volume

H0     = 1./(1 - Vi0 + ep_v.*V0vq + ep_d.*V0dq);  % constant in front

k1v     = 4.3.*nu0v.*E0v.*TE;
k1d     = 4.3.*nu0d.*E0d.*TE;
k2v     = ep_v.*r0v.*E0v.*TE;
k2d     = ep_v.*r0d.*E0d.*TE;
k3v     = 1 - ep_v;
k3d     = 1 - ep_d;
%%%%%
k1t     = 4.3.*nu0p.*E0p.*TE;
k1csf   = 0.04.*nu0p.^2.*E0p.^2.*TE;
k2p     = ep_p.*r0p.*E0p.*TE;
k3p     = 1 - ep_p;
wt      = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neuronal model parameters
%--------------------------------------------------------------------------


dt       = 0.005; % M.dt; % integration step (note must be smaller than 0.1 for more than 6 layers)
T        = 60/dt; %M.T;% time course lenght (sec/dt)
st_duration = 20;  % sec
% 
% %lam      = 0.3;                 % how slow is the inhibitory activity
% mu       = 2;                  % inhibitory to excitatory
% 
% A        = -3*eye(K);           % self-inhibition of excitatory
%A        = M.A;
 U        = zeros(T,K);
 U(5/dt:(5+st_duration)/dt,:) = 1.0;
Um   = zeros(T,1);
Um((5+st_duration-6)/dt:end,1)  = 1;
Um((5+st_duration)/dt+1:end,2)  = 1;

% Um((5+st_duration-6)/dt:(5+st_duration-5.5)/dt,1)  = 1;
% Um((5+st_duration)/dt+1:end,2)  = 1;
% x = zeros(2,1);
% Bsig = 0;
% Bmu = -1.5;
% Bin = 0;
% y_up  = x;
% for t = 1:T
% y_up(1) = y_up(1) + dt*(-(1.5+Bsig*Um(t,2))*x(1) - (1.5+Bmu*Um(t,2))*x(2) + U(t,1));
% y_up(2) = y_up(2) + dt*(0.1*(-x(2) + (1+Bin*Um(t,1))*x(1)));
% x = y_up;
% Y_up(t,:) = y_up; 
% end
% x = zeros(2,1);
% Bsig = -1*ones(1,1);
% Bmu = -1.5*ones(1,1);
% Bin = -20*ones(1,1);
% 
% 
% y_down  = x;
% for t = 1:T
% y_down(1) = y_down(1) + dt*(-(1.5+Bsig*Um(t,2))*x(1) - (1.5+Bmu*Um(t,2))*x(2) + U(t,1));
% y_down(2) = y_down(2) + dt*(0.1*(-x(2) + (1+Bin*Um(t,1))*x(1)));
% x = y_down;
% Y_down(t,:) = y_down; 
% end 
% figure, plot([Y_down(:,1),Y_up(:,1)]);
% %U        = M.U;
%  U        = zeros(T,K);
%  U(10/dt:(10+st_duration)/dt,:) = 1.0;
% x = zeros(K,2);


% y   = x;
% for t = 1:T
% y(:,1) = y(:,1) + dt*(-(1.5+Bsig*Um(t,2)).*x(:,1) - (1.5+Bmu*Um(t,2)).*x(:,2) + U(t,:)');
% y(:,2) = y(:,2) + dt*(0.1*(-x(:,2) + (1+Bin*Um(t,1)).*x(:,1)));
% x = y;
% Y(t,:) = y(:); 
% end 
%Um       = M.Um;
% NVC parameters:
%--------------------------------------------------------------------------
c1       = 0.6;
c2       = 1.5;
c3       = 2;
M.comp   = 5;
M.mean   = [0.1 0.25 0.5 0.75 0.9];
M.var    = [0.025 0.025 0.025 0.025 0.025];


M.amp    = [4. 7 3.5 7 4.]*0.1;
l        = linspace(0,1,K);
for i = 1:M.comp
    M.profile(:,i) = 1./sqrt(2*pi*M.var(i)).*exp(-(l-M.mean(i)).^2./(2.*M.var(i)));
    M.profile(:,i) = M.profile(:,i)./max(M.profile(:,i)).*M.amp(i);
end
C        = {diag(sum(M.profile,2))}; %


%==========================================================================
% SIMULATION SCENARIO: SIM 2 (Bump: (a) in Macrovasulature V0v (b) Neuronal (c) V0v bumb absolute F constant)
%==========================================================================



time_ax  = linspace(0,(T*dt),T)-5;
layer_ax = linspace(0,100,2*K+1);
layer_ax = layer_ax(2:2:end);



%cm    = jet(length(SimV0));
%scale = linspace(0.2,1,length(SimV0dr));

onset       = find(time_ax==min(abs(time_ax)));
offset      = find(abs(time_ax-st_duration)==min(abs(time_ax-st_duration)))-50;
l           = 0;



%BOLD_venules = [];
%BOLD_drein   = [];
%
% BOLD_EV      = [];
% BOLD_IV      = [];
% BOLD_volume  = [];
% BOLD2 = [];
% CMRO2 = [];
% FLOW  = [];
% dQd   =[];
% Vd    = [];
% Vv    = [];
% Va    = [];
% Qd = [];
% Qv = [];
% Qa = [];



X  = zeros(K,10);
Xp = zeros(1,2);
yp = Xp;
Yp = [];

y = X;
Y = [];

f_d = ones(K,1);
dv_d = ones(K,1);
dHb_d = ones(K,1);



for t = 1:T
    
    X(:,[4:10]) = exp(X(:,[4:10]));
    Xp         = exp(Xp);
    
    
    if t>(offset+50),
        
        mu = 0.5;
        lam      = 0.1;
    end
%     
    
    %--------------------------------------------------------------------------
    % Neuronal (excitatiry & inhibitory)
   % y(:,1)   = y(:,1) + dt*(A*X(:,1) - mu.*X(:,2) + C{1}*U(t,:)');
    y(:,1) = y(:,1) + dt*(-(1.5+Bsig*Um(t,1)).*X(:,1) - (1.5+Bmu*Um(t,2)).*X(:,2) + C{1}*U(t,:)');
    y(:,2) = y(:,2) + dt*(0.1*(-X(:,2) + (1+Bin*Um(t,1)).*X(:,1)));
    %        y(:,1)   = y(:,1) + dt*(A*X(:,1) - mu.*X(:,2) + C(:,siminp)*U(siminp,t));
   % y(:,2)   = y(:,2) + dt*(lam.*(-X(:,2) +  X(:,1)));
    %--------------------------------------------------------------------------
    % vasoactive signal:
    y(:,3)   = y(:,3) + dt*(X(:,1) - c1.*(X(:,3)));
    %--------------------------------------------------------------------------
    % Arterial (capilary) inflow:
    df_a     = c2.*X(:,3) - c3.*(X(:,4)-1);
    y(:,4)   = y(:,4) + dt*(df_a./X(:,4));
    %--------------------------------------------------------------------------
    
    % arterioles
   % f_a     = X(:,9).^(1./aa);
    % f_a     = (V0aq.*X(:,9).^(1./aa) + F0vq.*tau_a.*X(:,4))./(V0aq+F0aq.*tau_a);
  %  dv_a    = F0aq.*(X(:,4) - f_a)./V0aq;
  %  na      = 20;
  %  ma      =  (X(:,4)+na-1)./na;
   % CM2a(t,:) = ma;
  %  dHb_a   = F0aq.*(X(:,4) - f_a.*X(:,10)./X(:,9))./V0aq;
    
  %  y(:,9)  = y(:,9) + dt*(dv_a./X(:,9));
   % y(:,10)  = y(:,10)+ dt*(dHb_a./X(:,10));
    
    
    % VENULES COMPARTMENT:
    %--------------------------------------------------------------------------
    % blood outflow from venules compartment
    if sum(alpha_v)>0
        f_v     = (V0vq.*X(:,5).^(1./alpha_v) + F0vq.*tau_v.*X(:,4))./(V0vq+F0vq.*tau_v);
    else
        f_v     = f_a;
    end

    dv_v      = F0vq.*(X(:,4) - f_v)./V0vq;
    
    % derivative blood volume venunles:
    %--------------------------------------------------------------------------
    % rCMRO2:
    m        = (X(:,4) + n-1)./n;
    %
    CMRO2(t,:) = m;
    %-------------------------------------------------------------------------
    % derivative deoxyhemoglobin content draining vein:
    dHb_v      = F0vq.*(m - f_v.*X(:,6)./X(:,5))./V0vq;
    %      dHb_v      = (f_a.*F0aq.*X(:,10)./X(:,9) - f_v.*F0vq.*X(:,6)./X(:,5))./V0vq + ((X(:,10)./X(:,9)+X(:,6)./X(:,5))/2 - EF);
    %  f_v(i).*F0vq(i).*X(i,6)./X(i,5)
    %       net(t,:) = ((X(:,10)./X(:,9)+X(:,6)./X(:,5))/2 - EF);
    % blood volume update venules:
    y(:,5)     = y(:,5) + dt*(dv_v./X(:,5));
    % deoxyhemoglobin content update venules:
    y(:,6)     = y(:,6) + dt*(dHb_v./X(:,6));
    
    % blood outflow from draining vein compartment (deepest layer):
    if alpha_d(end)>0
        f_d(end)      = (V0dq(end).*X(end,7).^(1./alpha_d(end)) + tau_d(end).*f_v(end).*F0vq(end))./(V0dq(end)+F0dq(end).*tau_d(end));
    else
        f_d(end)      = f_v(end)*F0vq(end)./F0dq(end);
    end
    dv_d(end)     = (f_v(end).*F0vq(end) - f_d(end).*F0dq(end))./V0dq(end);
    dHb_d(end)    = (f_v(end).*F0vq(end).*X(end,6)./X(end,5) - f_d(end).*F0dq(end).*X(end,8)./X(end,7))./V0dq(end);
    
    % derivative blood volume draining vein (deepest layer):
    for i = K-1:-1:1,
        if alpha_d(i)>0
            f_d(i)     = (V0dq(i).*X(i,7).^(1./alpha_d(i)) + tau_d(i).*(f_v(i).*F0vq(i)+f_d(i+1).*F0dq(i+1)))./(V0dq(i)+F0dq(i).*tau_d(i));
        else
            f_d(i)     = f_v(i)*F0vq(i)./F0dq(i)+f_d(i+1)*F0dq(i+1)./F0dq(i);
        end
        
        dv_d(i)    = (f_v(i).*F0vq(i) + f_d(i+1).*F0dq(i+1) - f_d(i).*F0dq(i))./V0dq(i);
        dHb_d(i)   = (f_v(i).*F0vq(i).*X(i,6)./X(i,5) + f_d(i+1).*F0dq(i+1).*X(i+1,8)./X(i+1,7) - f_d(i).*F0dq(i).*X(i,8)./X(i,7))./V0dq(i);

    end;
    %Fout_d(t,:) = f_d;
    
    % blood volume update draining vein (deepest layer):
    y(:,7)  = y(:,7) + dt*(dv_d./X(:,7));
    y(:,8)  = y(:,8) + dt*(dHb_d./X(:,8));
    
    % pial vein:
    
    if alpha_p>0
        f_p     = (V0pq.*Xp(1).^(1./alpha_p) + F0pq.*tau_p.*f_d(1))./(V0pq+F0pq.*tau_p);
    else
        f_p     = f_d(1);
    end;
    %Fout_p(t,:) = f_p;
    dv_p  = (f_d(1).*F0dq(1) - f_p.*F0pq)./V0pq;
    dHb_p = (f_d(1).*F0dq(1).*X(1,8)./X(1,7) - f_p.*F0pq.*Xp(2)./Xp(1))./V0pq;
    
    
    yp(:,1)  = yp(:,1) + dt*(dv_p./Xp(1));
    yp(:,2)  = yp(:,2) + dt*(dHb_p./Xp(2));
    
    %                                               tau_v((X(:,4)-f_a)>=0) = tau_va((X(:,4)-f_a)>=0);
    %                                               tau_v((X(:,4)-f_a)<0)  = tau_vd((X(:,4)-f_a)<0);
%     if t<=(offset+50)
%         tau_v     = tau_v_in;
%         tau_d     = tau_d_in;
%         tau_p     = tau_p_in;
%     else
%         tau_v     = tau_v_de;
%         tau_d     = tau_d_de;
%         tau_p     = tau_p_de;
%     end
     tau_v     = tau_v_in;
    tau_d     = tau_d_in;
    tau_p     = tau_p_in;
 
    % check for deflation (negative derivative)
    tau_v(dv_v<0)  = tau_v_de(dv_v<0);
    tau_d(dv_d<0)  = tau_d_de(dv_d<0);
    tau_p(dv_p<0)  = tau_p_de(dv_p<0);
 
    
   % test(t,:)  = f_v - f_d;
   % test2(t,:) = f_a - f_v;
    X        = y;
    Xp       = yp;
   % Y(t,:)   = vec(y');
   % Yp(t,:)  = vec(yp');
   % Flow:
    NEU(t,:)   = y(:,1); 
    CBF(t,:)   = exp(y(:,4));
   % arterioles
   % v_a  = exp(y(:,9));
   % q_a  = exp(y(:,10));
    
    % venules:
    v_v  = exp(y(:,5));
    q_v  = exp(y(:,6));
    % draining vein:
    v_d  = exp(y(:,7));
    q_d  = exp(y(:,8));
    % pail vein:
    v_p  = exp(yp(:,1));
    q_p  = exp(yp(:,2));
    
    Qd(t,:) = q_d;
    Qv(t,:) = q_v;
    Qp(t,:) = q_p;

    Vd(t,:) = v_d;
    Vv(t,:) = v_v;
    Vp(t,:) = v_p;
%     Vi   = v_v.*V0v + v_d.*V0d;
    
    
    BOLD(t,:) = H0.*((1-Vi0).*( k1v.*V0vq.*(1-q_v)      +k1d.*V0dq.*(1-q_d))+...
                               +k2v.*V0vq.*(1-q_v./v_v) +k2d.*V0dq.*(1-q_d./v_d)+...
                               +k3v.*V0vq.*(1-v_v)      +k3d.*V0dq.*(1-v_d)).*100;
    
    
    %
  %  BOLD_arterioles(t,:) = (1./(1 - V0a + ep_a.*V0a)).*((1-V0a).*(k1a.*V0a.*(1-q_a))+k2a.*V0a.*(1-q_a./v_a)+k3a.*V0a.*(1-v_a)).*100;
    
    BOLD_venules(t,:) = (1./(1 - V0vq + ep_v.*V0vq)).*((1-V0vq).*(k1v.*V0vq.*(1-q_v))+k2v.*V0vq.*(1-q_v./v_v)+k3v.*V0vq.*(1-v_v)).*100;
    
    BOLD_drain(t,:)   = (1./(1 - V0dq + ep_d.*V0dq)).*((1-V0dq).*(k1d.*V0dq.*(1-q_d))+k2d.*V0dq.*(1-q_d./v_d)+k3d.*V0dq.*(1-v_d)).*100;
    BOLD_pail_CSF(t,:)    = (V0pq./(1 - V0pq + ep_p.*V0pq)).*((1-V0pq).*(wt.*k1t.*(1-q_p)+(1-wt).*k1csf.*(1-q_p.^2./v_p))+k2p.*(1-q_p./v_p)+k3p.*(1-v_p)).*100;
 % BOLD_pail_CSF(t,:)    = (V0pq./(1 - V0pq + ep_p.*V0pq)).*((1-V0pq).*(k1p.*V0pq.*(1-q_p))+k2p.*V0pq.*(1-q_p./v_p)+k3p.*V0pq.*(1-v_p)).*100;
     pail(t) = k1t.*(q_p-1);
     CSF(t) = k1csf.*(q_p.^2./v_p-1);
    %                                               %
  %  BOLD_EV(t,:)      = H0.*((1-Vi0).*(k1a.*V0a.*(1-q_a)+k1v.*V0v.*(1-q_v)+k1d.*V0d.*(1-q_d))).*100;
   % BOLD_IV(t,:)      = H0.*(k2a.*V0a.*(1-q_a./v_a)+k2v.*V0v.*(1-q_v./v_v)+k2d.*V0d.*(1-q_d./v_d)).*100;
   % BOLD_volume(t,:)  = H0.*(k3a.*V0a.*(1-v_a)+k3v.*V0v.*(1-v_v)+k3d.*V0d.*(1-v_d)).*100;
    
    
end



% figure,  plot([BOLD(:,1),BOLD_pail_CSF])
% 
% figure, mesh(NEU)
% figure, mesh(BOLD)