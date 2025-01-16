function est_kernel = BOLD_estimate_laminar_PSF(N,K,MAT_ind,depth_map_HR,depth_map,ind_vox_sel)
%load Data_laminar_noise;
% INPUT:
% N       - number of neuronal depths
% K       - number of BOLD depths
% Mat_ind - transformation from high-res to lower voxel resultion of
%           functional data
% depth_map_HR - high-resultion map of normalized corical depth [anatomy voxel space]
% depth_map    - lower-resultion map of normalized corical depth [functional voxel space]
% ind_vox_sel  - indices of selected voxels (e.g. indices of after specific ROI selection)
% 
% OUTPUT:
% est_kernel - estimated PSF for K-number of BOLD depths



m      = 0.5 + 0.06*randn(10,1);
nsig   = (1/(4*(N+1)))^2;


% Generate high res profiles in space for given curvature
IF = [];
for i = 1:length(m)
   IF(:,:,i) = 1./sqrt(2*pi*nsig).*exp(-(depth_map_HR-m(i)).^2./(2.*nsig));
end


IF = reshape(IF,[size(IF,1)*size(IF,2),size(IF,3)]);

[EV_vox_s,ind0] = unique(depth_map_HR(:));

% downsample
IFr = zeros(length(unique(MAT_ind(:))),size(IF,2));
for j = 1:size(IF,2)
    for i = 1:length(unique(MAT_ind(:)))
        IFr(i,j)  = nanmean(IF(MAT_ind==i,j));
    end
end

[dist, ind] = sort(depth_map(ind_vox_sel));

yy = IFr(ind_vox_sel,:);
yy = IFr(ind_vox_sel,:);
yy = yy(ind,:);

IFm = interp1(EV_vox_s,IF(ind0,:),dist);
IFm(isnan(IFm)) = 0;

params(1) = 1;
params(2) = 0.1;

inp = IFm;
est = fminsearch(@(p)ConvKernel(p,yy,inp,dist),params);


[err, kernel yp sp] = ConvKernel(est,yy,inp,dist);

layer_axr = linspace(0,1,2*K+1);
layer_axr = layer_axr(2:2:end);


est_kernel = interp1(sp,kernel,layer_axr,'pchip','extrap');

est_kernel = est_kernel./sum(est_kernel);


return
