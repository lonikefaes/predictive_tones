function V = surf2vox(FV,range,interp_factor)



FV.vertices = interp_factor/2*(FV.vertices-repmat([range.r2.mx(1),range.r2.my(1),range.r2.mz(1)],size(FV.vertices,1),1));

nx = length(range.r5.mx);
ny = length(range.r5.my);
nz = length(range.r5.mz);

V = VOXELISE(1:nx,1:ny,1:nz,FV);


return
% V0 = zeros(nx,ny,nz);
% V0(4:end-3,4:end-3,4:end-3) = 1;

[xx,yy,zz] = ndgrid(-1:1);
nhood      = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1;

% 2step procedure
%V3 = logical(V3 + ~V0).*~logical(V1+V2);
V1 = imerode(imdilate(V10,nhood),nhood);
%V2 = imerode(imdilate(V20,nhood),nhood);
V2 = V20;
V3 = imerode(imdilate(V30,nhood),nhood);
V1p = V1(:,:,interp_factor+1:end-interp_factor);
V2p = V2(:,:,interp_factor+1:end-interp_factor);
V3p = V3(:,:,interp_factor+1:end-interp_factor);
% for i = 1:1
%     
%     V3p = imdilate(V3p,nhood).*~V2p.*~V1p;
% %     if i<3
% %         V1p = imdilate(V1p,nhood).*~V2p;
% %         V2p = imdilate(V2p,nhood).*~V1p.*~V3p;
% %     end
% end
V1(:,:,interp_factor+1:end-interp_factor) = V1p;
V2(:,:,interp_factor+1:end-interp_factor) = V2p;
V3(:,:,interp_factor+1:end-interp_factor) = V3p;
BWV1 = bwdist(V1);
BWV2 = bwdist(V2);
BWV3 = bwdist(V3);%.*~V2.*~V1;
Mask1 = (BWV1-BWV2)>0;
Mask2 = (BWV2-BWV3)>0;% - Mask1; % GM
Mask3 = (BWV1-BWV3)>0;
Mask1(Mask2|Mask3) = 1;
Mask  = Mask1 + 1*Mask2; % WM = 2, GM = 1, CSF = 0; 
%Mask(:,:,interp_factor+1:end-interp_factor) = 2*V3p + V2p;


figure(1);
for k = 1:nz
   subplot(121); imagesc(Mask(:,:,k)); colorbar; axis image; drawnow; title(num2str(k));
   subplot(122); imagesc(V10(:,:,k)+V20(:,:,k)*2+V30(:,:,k)*3); colorbar; axis image; drawnow; title(num2str(k));
   pause;
end

pause;
Mask3 = Mask==2; % WM
Mask2 = Mask==1; % GM
Mask1 = Mask==0; % CSF




[X0 Y0 Z0]  = meshgrid(0:ny0-1,0:nx0-1,0:nz0-1);
[Xr Yr Zr]  = meshgrid(linspace(0,ny0-1,ny0*interp_factor/2),...
                       linspace(0,nx0-1,nx0*interp_factor/2),...
                       linspace(0,nz0-1,nz0*interp_factor/2));
                     
Act  = Act(mx2,my2,mz2);
ActR = interp3(X0,Y0,Z0,Act,Xr,Yr,Zr,'nearest');
Ant  = Ant(mx2,my2,mz2);
AntR = interp3(X0,Y0,Z0,Ant,Xr,Yr,Zr,'nearest');
UR   = interp3(X0,Y0,Z0,U,Xr,Yr,Zr,'linear');
Vm   = Vm(mx0,my0,mz0);