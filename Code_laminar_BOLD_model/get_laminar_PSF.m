function kernel = get_laminar_PSF(N,K,MAT_ind,EV_vox,EV_vox_d,ind_vox,sel_sel,Hs)
%load Data_laminar_noise;



depths = 30;
m      = 0.5 + 0.06*randn(10,1);
nsig   = (1/(4*(N+1)))^2;


% Generate high res profiles in space for given curvature
IF = [];
for i = 1:length(m)
   IF(:,:,i) = 1./sqrt(2*pi*nsig).*exp(-(EV_vox-m(i)).^2./(2.*nsig));
end


IF = reshape(IF,[size(IF,1)*size(IF,2),size(IF,3)]);

[EV_vox_s,ind0] = unique(EV_vox(:));


% downsample
interp_factor = 6;
IFr = zeros(length(unique(MAT_ind(:))),size(IF,2));
for j = 1:size(IF,2)
    for i = 1:length(unique(MAT_ind(:)))
        IFr(i,j)  = nanmean(IF(MAT_ind==i,j));
    end
end

[dist ind] = sort(EV_vox_d(ind_vox(sel_sel)));
[dist0v ind0v] = sort(EV_vox_d);
yy = IFr(ind_vox(sel_sel),:);
yy = IFr(ind_vox(sel_sel),:);
yy = yy(ind,:);
%sp = [0:0.01:1]';
IFm = interp1(EV_vox_s,IF(ind0,:),dist);
IFm(isnan(IFm)) = 0;
%yy = interp1(dist,yy,sp,'pchip','extrap');
figure(3),
plot(dist,yy,'.',EV_vox_s,IF(ind0,:),'MarkerSize',20); axis square; xlabel('Cortical depth'); ylabel('(a.u)')
xlim([0 1]), ylim([0 8])
hold on;
% Design convlution kernel to explain difference between n2k an yd!!!
% nonlinear least squares
params(1) = 1;
params(2) = 0.1;
%params(3) = 1;
inp = IFm;
est = fminsearch(@(p)ConvKernel(p,yy,inp,dist),params);

%modelFun =  @(p,x) exp(-0.5.*((x(:)-0.5)./p(2)).^2).^p(3);



[err, kernel yp sp] = ConvKernel(est,yy,inp,dist);

plot(dist,yp);
ax=gca;ax.XTick = [0 0.5 1]; ax.XTickLabel = {'WM','Mid-GM','CSF'}; xlim([0 1]); hold off
kernel2  = interp1(sp,kernel,dist,'pchip','extrap');
H = [];
% layer_ax = linspace(0,1,2*(K)+1);
% layer_ax = layer_ax(2:2:end);
layer_axr = linspace(0,1,2*K+1);
layer_axr = layer_axr(2:2:end);

layer_axr_ex = [layer_axr(1)-(layer_axr(2)-layer_axr(1)),layer_axr,layer_axr(end)+(layer_axr(2)-layer_axr(1))];
for i = 1:length(layer_axr_ex),
   H(:,i) = (exp(-(length(layer_axr_ex).*(dist-layer_axr_ex(i))).^2)).^3;
end
H  = H./repmat(sum(H,2),1,length(layer_axr_ex));
H  = H(:,2:end-1);
Hs = H./repmat(sum(H,1),size(H,1),1);



kernel2 = Hs'*kernel2;
kernel2 = kernel2./sum(kernel2);

kernel3  = interp1(sp,kernel,layer_axr,'pchip','extrap');

% m      = linspace(0,1,2*N+1);
% m      = repmat(m(2:2:end),K,1);
% nsig   = (1/(4*(N+1)))^2;
% sp     = linspace(0,1,2*K+1);
% sp     = repmat(sp(2:2:end)',1,N);
% 
% n2k    = 1./sqrt(2*pi*nsig).*exp(-(sp-m).^2./(2.*nsig));
% 
% figure(5), plot(layer_axr,[n2k,yd]);
% 
% % Design convlution kernel to explain difference between n2k an yd!!!
% % nonlinear fitting
% 
% 
% options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxFunEvals',10000,'TolFun',1e-12,'MaxIter',2000);
% 
% % nonlinear least squares
% params(1) = 0;
% params(2) = 1;
% params(3) = 0;
% est = lsqnonlin(@ConvKernel,params,[],[],options,yd,n2k,N,K); 
% 
% [err, kernel yp] = ConvKernel(est,yd,n2k,N,K);
% figure(6), plot(layer_axr,[sum(yp,2),sum(yd,2)]);
kernel3 = kernel3./sum(kernel3);
figure(7), plot(layer_axr,kernel2,layer_axr,kernel3,sp,kernel*max(kernel3)); axis square;
ax = gca; ax.XTick = [0 0.5 1]; ax.XTickLabel = {'WM','Mid-GM','CSF'}; xlim([0 1]); hold off
kernel = kernel3;
return
% 
% layer_axr = linspace(0,1,2*(K)+1);
% layer_axr = layer_axr(2:2:end);
% 
% H = [];
% IF = [];
% layer_axr_ex = [layer_axr(1)-(layer_axr(2)-layer_axr(1)),layer_axr,layer_axr(end)+(layer_axr(2)-layer_axr(1))];
% for i = 1:length(layer_axr_ex),
%    H(:,i) = (exp(-(length(layer_axr_ex).*(dist-layer_axr_ex(i))).^2)).^3;
%    IF(:,:,i) = exp(-(length(layer_axr_ex).*(EV_vox-layer_axr_ex(i))).^2);
%    n2k    = 1./sqrt(2*pi*nsig).*exp(-(sp-m).^2./(2.*nsig));
% end
% IF = IF(:,:,2:end-1);
% IF = reshape(IF,[size(IF,1)*size(IF,2),size(IF,3)]);
% interp_factor = 6;
% IFr = zeros(length(unique(MAT_ind(:))),size(IF,2));
% for j = 1:size(IF,2)
%     for i = 1:length(unique(MAT_ind(:)))
%         IFr(i,j)  = nansum(IF(MAT_ind==i,j))./(interp_factor.^2);
%     end
% end
% 
% figure, plot(EV_vox(:),IF(:,1),'.')
% [Ds,inds] = sort(EV_vox(:));
% IFs = IF(inds,5);
% [Ds,inds]= unique(Ds);
% IFs = IFs(inds);
% IFF = interp1(Ds,IFs,[0:0.01:1]');
% figure, plot(IFF)
