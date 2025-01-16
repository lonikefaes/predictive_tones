function [div3,Gx,Gy,Gz] = graddiv(potential,mask)

[nx, ny, nz]   = size(potential);
[gx, gy, gz]   = gradnan(potential);


gn    = sqrt(gx.^2+gy.^2+gz.^2);
I     = find(gn<eps);
gn(I) = 1;
gx    = gx./gn;
gy    = gy./gn;
gz    = gz./gn;

%fwhm    = [fwhm 30];

% Number of basis functions for x, y & z
%-----------------------------------------------------------------------
% 
% potential = zeros(75,40);
% 
% potential(30:55,:)  = repmat(linspace(0.01,0.99,26)',1,40);
% potential(56:end,:) = 1;
% 
% [nx ny nz] = size(potential);
cutoff = 6;
tmp    = sqrt(sum(eye(3).^2));
k      = max(round(([nx,ny,nz].*tmp)/cutoff),[1 1 1]);

% Scaling is to improve stability.
%-----------------------------------------------------------------------
stabilise = 8;
basX = spm_dctmtx(nx,k(1))*stabilise;
basY = spm_dctmtx(ny,k(2))*stabilise;
basZ = spm_dctmtx(nz,k(3))*stabilise;

dbasX = spm_dctmtx(nx,k(1),'diff')*stabilise;
dbasY = spm_dctmtx(ny,k(2),'diff')*stabilise;
dbasZ = spm_dctmtx(nz,k(3),'diff')*stabilise;

% d2basX = spm_dctmtx(nx,k(1),'diff2')*stabilise;
% d2basY = spm_dctmtx(ny,k(2),'diff2')*stabilise;
% d2basZ = spm_dctmtx(nz,k(3),'diff2')*stabilise;

BAS = kron(basZ,kron(basY,basX));
%BAS = kron(basY,basX);
%dBAS = kron(dbasZ,kron(dbasY,dbasX));
%dBAS = kron(dbasY,dbasX);

potential = potential./max(potential(:));
% B = BAS\potential(:);
% P = zeros(nx,ny,nz);
% P(:) = BAS*B; 

% 
% % X-gradient
% 
% for i = 1:ny
%     for j = 1:nz
%       betaX     = basX\P(:,i,j);
%       Px(:,i,j) = basX*betaX; 
%       dPx(:,i,j) = dbasX*betaX;
%       d2Px(:,i,j) = d2basX*betaX;
%     end
% end
% for i = 1:ny
%     for j = 1:nz
%       betaX     = basX\P(:,i,j);
%       Px(:,i,j) = basX*betaX; 
%       dPx(:,i,j) = dbasX*betaX;
%       d2Px(:,i,j) = d2basX*betaX;
%     end
% end
% % Y gradient
% for i = 1:nx
%     for j = 1:nz
%       betaY     = basY\P(i,:,j)';
%       Py(i,:,j) = basY*betaY; 
%       dPy(i,:,j) = dbasY*betaY;
%       d2Py(i,:,j) = d2basY*betaY;
%     end
% end
% % Z gradient
% for i = 1:nx
%     for j = 1:ny
%       betaZ     = basZ\squeeze(P(i,j,:));
%       Pz(i,j,:) = basZ*betaZ; 
%       dPz(i,j,:) = dbasZ*betaZ;
%       d2Pz(i,j,:) = d2basZ*betaZ;
%     end
% end
% figure,
% k = 10;
% Px(mask~=1) = NaN;
% dPx(mask~=1) = NaN;
% d2Px(mask~=1) = NaN;
% dPy(mask~=1) = NaN;
% d2Py(mask~=1) = NaN;
% dPz(mask~=1) = NaN;
% d2Pz(mask~=1) = NaN;
% subplot(231); imagesc(d2Px(:,:,k)); colorbar; subplot(232); imagesc(d2Py(:,:,k)); colorbar; subplot(233); imagesc(d2Pz(:,:,k)); colorbar;
% subplot(234); imagesc(dPx(:,:,k)); colorbar; subplot(235); imagesc(dPy(:,:,k)); colorbar; subplot(236); imagesc(dPz(:,:,k)); colorbar;
% 
% divP = d2Px + d2Py + d2Pz;
% 
% figure,
% imagesc(divP(:,:,15)); colorbar

%Bx = vb_glm(gx(mask==1),BAS(vec(mask==1),:),0);
Bx = BAS(vec(mask==1),:)\gx(mask==1);
By = BAS(vec(mask==1),:)\gy(mask==1);
%By = vb_glm(gy(mask==1),BAS(vec(mask==1),:),0);
Bz = BAS(vec(mask==1),:)\gz(mask==1);
%Bz = vb_glm(gz(mask==1),BAS(vec(mask==1),:),0);

% 
% G  = zeros(nx,ny);
% dG = zeros(nx,ny);
% beta  = BAS\potential(:);
% G(:)  = BAS*beta;
% dG(:) = dBAS*beta;
% dG(:) = BAS(:,2:8:end)*beta(2:8:end);
% figure,
% subplot(121);
% imagesc(G), colorbar;
% subplot(122), imagesc(dG); colorbar;


%G = zeros(nx,ny,nz);
%GGx = zeros(nx,ny,nz);

Gx = zeros(nx,ny,nz);
Gy = zeros(nx,ny,nz);
Gz = zeros(nx,ny,nz);
%G(:)        = BAS*B;
%GGx(:)      = BAS(:,1:5:end)*B(1:5:end);
Gx(mask==1) = BAS(mask==1,:)*Bx;
Gy(mask==1) = BAS(mask==1,:)*By;
Gz(mask==1) = BAS(mask==1,:)*Bz;

% Gx(:) = BAS(:,:)*Bx;
% Gy(:) = BAS(:,:)*By;
% Gz(:) = BAS(:,:)*Bz;
% figure, 
% for k = 1:nz
% subplot(221); imagesc(potential(:,:,k).*(mask(:,:,k)==1)); colorbar; axis image; subplot(222); imagesc(G(:,:,k).*(mask(:,:,k)==1)); impixelinfo; colorbar; axis image; 
% subplot(223); imagesc(GGx(:,:,k).*(mask(:,:,k)==1)); colorbar; axis image; 
% pause
% end
% k = 10;
% figure, subplot(2,3,1), imagesc(Gx(:,:,k)); colorbar; axis image; subplot(2,3,4), imagesc(gx(:,:,k)); colorbar; axis image
%         subplot(2,3,2), imagesc(Gy(:,:,k)); colorbar; axis image; subplot(2,3,5), imagesc(gy(:,:,k)); colorbar; axis image
%         subplot(2,3,3), imagesc(Gz(:,:,k)); colorbar; axis image; subplot(2,3,6), imagesc(gz(:,:,k)); colorbar; axis image
%         



% gx(mask~=1) =  NaN;
% gy(mask~=1) =  NaN;
% gz(mask~=1) =  NaN;

Gn    = sqrt(Gx.^2+Gy.^2+Gz.^2);
I     = find(Gn<eps);
Gn(I) = 1;
Gx    = Gx./Gn;
Gy    = Gy./Gn;
Gz    = Gz./Gn;

% calculate divergence;
[px, ~, ~] = gradnan(Gx);
[~, qy, ~] = gradnan(Gy);
[~, ~, rz] = gradnan(Gz);
div1       = px+qy+rz;
% div = divergence(gx,gy,gz);
% div1(isnan(div1)) = 0;

Gx(mask~=1) =  NaN;
Gy(mask~=1) =  NaN;
Gz(mask~=1) =  NaN;

[px, ~, ~] = gradnan(Gx);
[~, qy, ~] = gradnan(Gy);
[~, ~, rz] = gradnan(Gz);
div2       = px+qy+rz;
% figure,
% div1(mask~=1) = NaN;
% subplot(131);
% imagesc(div1(:,:,k)); colorbar
% subplot(132);
% imagesc(div(:,:,k)); colorbar
% subplot(133);
% imagesc(div2(:,:,k)); colorbar

%D = vb_glm(div2(mask==1 & ~isnan(div2)),BAS(mask==1 & ~isnan(div2),:),0);
D = BAS(mask==1 & ~isnan(div2),:)\vec(div2(mask==1 & ~isnan(div2)));
div3 = zeros(nx,ny,nz);
div3(mask==1) = BAS(mask==1,:)*D;
%div3(mask~=1) = ;
%figure, subplot(121); imagesc(div3(:,:,15)); colorbar; subplot(122); imagesc(div2(:,:,15)); colorbar;
%return