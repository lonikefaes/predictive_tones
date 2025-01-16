function [Error kernel yp sp] = ConvKernel(par,y,inp,dist)



% m      = 0.5;
% nsig0   = (1/(4*(N+1)));%^2;
% nsig   = nsig0*exp(par(1));
% if mod(K,2)==0
%     sp     = linspace(0,1,K-1+2);
% else
%     sp     = linspace(0,1,K+2);
% end
sp     = linspace(0,1,1000)';%dist;


%A      = par(1);

%e      = 1*exp(par(3));

%kernel = A./sqrt(2*pi*nsig).*(exp(-(sp-m).^2./(2.*nsig))).^e;
%kernel = A*(exp(-(nsig.*(sp-m)).^2)).^e;

%kernel = gbellmf(sp,par(1:3));

modelFun =  @(p,x) 1./sqrt(2*pi*p(2)).*exp(-(x(:)-0.5).^2./(2*p(2)));%.^p(3);

 %  IF(:,:,i) = 1./sqrt(2*pi*nsig).*exp(-(EV_vox-m(i)).^2./(2.*nsig));
%p0 = [4 0.1 1];


kernel  = modelFun(par,sp);
kernel  = kernel./sum(kernel);
inp = interp1(dist,inp,sp,'pchip','extrap');
%y2 = interp1(dist,y,sp,'pchip','extrap');
%inp = inp./max(inp);

%[q r]= deconv(inp,yp);
yp     = conv2(inp,kernel(:),'same'); 
yp     = interp1(sp,yp,dist);
yp     = par(1)*yp;
%yp2     = conv2(inp,r./sum(r),'same'); 
Error  = norm(yp(:)-y(:));
kernel  = kernel./max(kernel); 
%kernel  = interp1(sp,kernel,dist,'pchip','extrap');
