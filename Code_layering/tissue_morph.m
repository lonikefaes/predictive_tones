function Mask_morph = tissue_morph(Vox,interp_factor,mask)

CSF0 = Vox.csf;
GM0  = Vox.gm;
WM0  = Vox.wm;

[nx, ny,nz] = size(CSF0);

[xx,yy,zz] = ndgrid(-1:1);
nhood      = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1;

% 2step procedure
%V3 = logical(V3 + ~V0).*~logical(V1+V2);
CSF = imerode(imdilate(CSF0,nhood),nhood);
%V2 = imerode(imdilate(V20,nhood),nhood);
GM  = GM0;
WM  = imerode(imdilate(WM0,nhood),nhood);
% V1p = V1(:,:,interp_factor+1:end-interp_factor);
% V2p = V2(:,:,interp_factor+1:end-interp_factor);
% V3p = V3(:,:,interp_factor+1:end-interp_factor);
% % for i = 1:1
% %     
% %     V3p = imdilate(V3p,nhood).*~V2p.*~V1p;
% % %     if i<3
% % %         V1p = imdilate(V1p,nhood).*~V2p;
% % %         V2p = imdilate(V2p,nhood).*~V1p.*~V3p;
% % %     end
% % end
% V1(:,:,interp_factor+1:end-interp_factor) = V1p;
% V2(:,:,interp_factor+1:end-interp_factor) = V2p;
% V3(:,:,interp_factor+1:end-interp_factor) = V3p;
BW_CSF = bwdist(CSF);
BW_GM = bwdist(GM);
BW_WM = bwdist(WM);%.*~V2.*~V1;
Mask0 = zeros(nx,ny,nz);
% Mask0([1 nx],:,:) = 1;
% Mask0(:,[1 ny],:) = 1;
% Mask0(:,:,[1 nz]) = 1;

Mask0 = logical(Mask0 + mask);

Mask1 = (BW_CSF-BW_GM)>0;
Mask2 = (BW_GM-BW_WM)>0;% - Mask1; % GM
Mask3 = (BW_CSF-BW_WM)>0;
Mask1(Mask2|Mask3) = 1;
Mask_morph  = (Mask1 + Mask2).*~Mask0 + Mask0*3; % WM = 2, GM = 1, CSF = 0; 




return
pause;
Mask3 = Mask==2; % WM
Mask2 = Mask==1; % GM
Mask1 = Mask==0; % CSF