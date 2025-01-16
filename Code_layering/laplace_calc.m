function [U_out, U_outS, Ux] = laplace_calc(Mask_morph)

[nx,ny,nz] = size(Mask_morph);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization for Laplace solution:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b1 = 10000; % CSF boundary 
b2 = 0.0001; % GW boundary

U        = zeros(nx,ny,nz);
U(Mask_morph==0) = b1;
U(Mask_morph==2) = b2;
U(Mask_morph==1) = -Inf;           % -INF is needed where is GM
U(Mask_morph==3) = NaN;
Ux                       = ones(nx+2,ny+2,nz+2)*NaN; % Box around with insulation
Ux(2:nx+1,2:ny+1,2:nz+1) = U;

% calculate Laplace distribution (takes a minute...)
[U_out] = laplace3d(Ux); % direct solution via matrix inversion
%[U_out] = laplace3djacobi(Ux,tol); % This is iteative solution if you have
%many voxels - like from entire brain (but might need some adjustment the tolerance threshold...I don't use this one)
% remove the outside box
U_out  = U_out(2:end-1,2:end-1,2:end-1);

% Do it again but this time with slightly dilated GM (needed for better gradients and divergence, to avoid error at the edges of the GM)
[xx,yy,zz] = ndgrid(-1:1);
nhood      = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1;

Mask_gm = imdilate(Mask_morph==1,nhood);
Mask_csf = logical(Mask_morph==0.*~Mask_gm);
Mask_wm = logical(Mask_morph==2.*~Mask_gm);

% same thing as above
U2           = zeros(nx,ny,nz);
U2(Mask_csf) = b1;
U2(Mask_wm)  = b2;
U2(Mask_gm)  = -Inf;
U2(Mask_morph==3) = NaN;
Ux2                       = ones(nx+2,ny+2,nz+2)*NaN;
Ux2(2:nx+1,2:ny+1,2:nz+1) = U2;


[U_outS] = laplace3d(Ux2);
U_outS   = U_outS(2:end-1,2:end-1,2:end-1);