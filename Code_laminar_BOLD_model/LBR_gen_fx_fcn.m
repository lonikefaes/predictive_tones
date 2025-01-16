function [ht_f, ht_dfdx] = LBR_gen_fx_fcn(N,K)

% N = 3;   % no. of neuronal depths
% K = 12;  % no. of hemodynamic depths

% Neuronal and NVC model specification
A  = sym('a',[N,N]);
MU = sym('mu',[N,1]);
LAM = sym('lam',[N,1]);
C1 = sym('c1',[N,1]);
C2 = sym('c2',[N,1]);
C3 = sym('c3',[N,1]);
Xn = sym('xn',[N,4]);
CU  = sym('cu',[N,1]);
%U  = sym('u',[N,1]);


fn = [A*Xn(:,1) - MU.*Xn(:,2) + CU;
      LAM.*(-Xn(:,2) + Xn(:,1));
      Xn(:,1) - C1.*Xn(:,3);
      (C2.*Xn(:,3) - C3.*(Xn(:,4)-1))./Xn(:,4)];

% Hemodynamic model specification  
Xk = sym('xk',[K,4]);
AL_v = sym('al_v',[K,1]);
AL_d = sym('al_d',[K,1]);
TT_v = sym('tt_v',[K,1]);
TT_d = sym('tt_d',[K,1]);
VE_v = sym('ve_v',[K,1]);
VE_d = sym('ve_d',[K,1]);
NR   = sym('nr',[K,1]);
F0v  = sym('f0v',[K,1]);
F0d  = sym('f0d',[K,1]);
N2K  = sym('n2k',[K,N]);
% loop across cortical depths
for k = K:-1:1;
    Xnk(k)     = N2K(k,:)*(Xn(:,4)-1) + 1;
    fv(k)      = (TT_v(k).*Xk(k,1).^(1./AL_v(k)) + VE_v(k).*Xnk(k))./(TT_v(k)+VE_v(k));
    dv_v(k,1)  = (Xnk(k) - fv(k))./(TT_v(k)*Xk(k,1));
    dHb_v(k,1) = ((Xnk(k) + NR(k)-1)./NR(k) - fv(k).*Xk(k,2)./Xk(k,1))./(TT_v(k)*Xk(k,2));
    if k == K
        fd(k)      = (TT_d(k).*Xk(k,3).^(1./AL_d(k)) + VE_d(k).*fv(k))./(TT_d(k)+VE_d(k));
        dv_d(k,1)  = (fv(k).*F0v(k)./F0d(k) - fd(k))./(TT_d(k)*Xk(k,3));
        dHb_d(k,1) = (fv(k).*Xk(k,2)./Xk(k,1) - fd(k).*Xk(k,4)./Xk(k,3))./(TT_d(k)*Xk(k,4));
    else
        fd(k)      = (TT_d(k).*Xk(k,3).^(1./AL_d(k)) + VE_d(k).*(fv(k).*F0v(k)/F0d(k)+fd(k+1).*F0d(k+1)/F0d(k)))./(TT_d(k)+VE_d(k));
        dv_d(k,1)  = (fv(k).*F0v(k)./F0d(k) + fd(k+1).*F0d(k+1)./F0d(k) - fd(k))./(TT_d(k)*Xk(k,3));
        dHb_d(k,1) = (fv(k).*F0v(k)./F0d(k).*Xk(k,2)./Xk(k,1) + fd(k+1).*F0d(k+1)./F0d(k).*Xk(k+1,4)./Xk(k+1,3) - fd(k).*Xk(k,4)./Xk(k,3))./(TT_d(k)*Xk(k,4));
    end
end

fk = [dv_v;dHb_v;dv_d;dHb_d];


f_s    = [fn;fk];
dfdx_s = jacobian([fn;fk],[Xn(:);Xk(:)]); 

% fx function handle
ht_f = matlabFunction(f_s,'File',['sym_fx_LBR_N',num2str(N),'K',num2str(K)],...
                          'Vars',[Xn(:);Xk(:);A(:);MU(:);CU(:);LAM(:);C1(:);C2(:);C3(:);...
                                  N2K(:);AL_v(:);AL_d(:);TT_v(:);TT_d(:);NR(:);VE_v(:);VE_d(:);F0v(:);F0d(:)]);

% Jacobian handle                                                  
ht_dfdx = matlabFunction(dfdx_s,'File',['sym_dfdx_LBR_N',num2str(N),'K',num2str(K)],...
                                'Vars',[Xn(:);Xk(:);A(:);MU(:);CU(:);LAM(:);C1(:);C2(:);C3(:);...
                                        N2K(:);AL_v(:);AL_d(:);TT_v(:);TT_d(:);NR(:);VE_v(:);VE_d(:);F0v(:);F0d(:)]);

return                                    
% Test:
dt = 0.01;
T   = 40/dt;
P.A = eye(N)*-0.5;
P.C = eye(N)/4;
P.U = zeros(N,T);
P.U(:,2/dt:5/dt) = 1;
P.MU = ones(N,1);
P.LAM = ones(N,1)*0.1;
P.C1 = ones(N,1)*0.6;
P.C2 = ones(N,1)*1.5;
P.C3 = ones(N,1)*0.6;

V0t    = 3.5;
s_v    = 0;
s_d    = 0.4;
w_v    = 0.5;
depths = linspace(0,100,2*K+1); % Normalized distance to the center of individual depths (in %)l    = depths(2:2:end);
l    = depths(2:2:end);
x_v  = 10+s_v*flipud(l(:));
x_d  = 10+s_d*flipud(l(:));
V0v  = x_v./sum(x_v)*w_v*V0t;
V0d  = x_d./sum(x_d)*(1-w_v)*V0t;
P.AL_v = ones(K,1)*0.3;
P.AL_d = ones(K,1)*0.15;
P.TT_v = ones(K,1)*1;
P.F0v  = V0v./P.TT_v;
P.F0d  = flipud(cumsum(flipud(P.F0v)));
P.TT_d = V0d./P.F0d;
P.VE_v = ones(K,1)*2;
P.VE_d = ones(K,1)*20;
P.NR   = ones(K,1)*3;
P.N2K  = kron(eye(N),ones(K/N,1));


x  = zeros(N*4+K*4,1);
tic
X = zeros(T,N*4+K*4);
xn  = zeros(N,4);
xk  = zeros(K,4);
for t = 1:size(P.U,2)
    xn(:)   = x(1:N*4);
    xk(:)   = x(N*4+1:end);
    xn(:,4) = exp(xn(:,4));
    xk      = exp(xk);
    IN       = num2cell([xn(:);xk(:);P.A(:);P.MU(:);P.C*P.U(:,t);P.LAM(:);P.C1(:);P.C2(:);P.C3(:);...
                         P.N2K(:);P.AL_v(:);P.AL_d(:);P.TT_v(:);P.TT_d(:);P.NR(:);P.VE_v(:);P.VE_d(:);P.F0v(:);P.F0d(:)]);
    fx      = ht_f(IN{:});
    dfdx    = ht_dfdx(IN{:});
    x       = x(:) + dt*fx(:) + 0.5*dt^2*(dfdx*fx(:));
    X(t,:)  = x(:);
end
toc;

figure, plot(X);
