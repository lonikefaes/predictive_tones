function [EDV2, EVV2] = layer_dist(U_out,U_outS,Mask_morph,Ux)

b1 = 10000; % CSF boundary 
b2 = 0.0001; % GW boundary

[nx,ny,nz] = size(U_out);
% calculate gradients on Laplace solution for dilated GM volume
% there are NaN's so avoid them using this gradnan function 
[gx, gy, gz] = gradnan(U_outS);


% normalize gradients:
gn    = sqrt(gx.^2+gy.^2+gz.^2);
I     = find(gn<eps);
gn(I) = 1;
gx    = gx./gn;
gy    = gy./gn;
gz    = gz./gn;


% calculate divergence;
[px, ~, ~] = gradnan(gx);
[~, qy, ~] = gradnan(gy);
[~, ~, rz] = gradnan(gz);
div1       = px+qy+rz;
% 


div1(isnan(div1)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of streamline caluculation (Equivolume&Equidistant representation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%U(:,:,:) = 0; % don't consider bottom and top slices...

% select voxels that has to be evaluated (within GM)
UU          =  zeros(nx,ny,nz);
UU(Mask_morph==1) = 1;
UU0         = ones(nx,ny,nz);
UU0(UU0==1) = [1:length(find(UU0==1))];
UU0         = UU.*UU0;
UU(UU==1)   = [1:length(find(UU==1))];
% and get their indices and subcoordinates
ind         = UU0(UU>0);
ind0        = ind;
[vX,vY,vZ]  = ind2sub([nx ny nz],ind);

VectorF                = cat(4,gy,gx,gz);
opts.FiberLengthMax    = 300;
opts.FiberLengthMin    = 3;
opts.DeviationAngleMax = 15; %(degrees)
opts.Step              = 0.2;
FiberLengthMax         = 300;
FiberLengthMin         = 3;
DeviationAngleMax      = 20; %(degrees)
Step                   = 0.2;
Roi                    = zeros(nx,ny,nz);
Roi(Mask_morph==1)     = 1;

EV = zeros(length(ind),1);
ED = zeros(length(ind),1);
TC = zeros(length(ind),1);
FB = cell(length(ind),1);

try
    cores = feature('numcores');
    matlabpool(cores);
end

disp('Starting parallel process calculation');
blocks = round(linspace(1,size(vX,1),100));
h0=tic;
for p = 1:length(blocks)-1,
    parfor j = blocks(p):blocks(p+1),
        xyz = [vX(j),vY(j),vZ(j)];
        [EV(j,:), ED(j,:), TC(j,:), FB{j}] = evaluatestreamline3D(xyz,U_out,VectorF,div1,Roi,b1,b2,opts);
    end
    h1=toc(h0);
    disp(['Finished: ', num2str(100*blocks(p+1)/size(vX,1)),'%, ', num2str(h1/60),' minutes elapsed...']);
end
 % EQUIVOLUME:
 EVV = zeros(nx,ny,nz);
 EVV(ind0) = EV(1:length(ind0)); 
 % EQUIDISTANT:
 EDV = zeros(nx,ny,nz);
 EDV(ind0) = ED(1:length(ind0));
 
 Thickness       = zeros(nx,ny,nz);
 Thickness(ind0) = TC(1:length(ind0)); 
  
%  figure,
%  
%  for k = 1:nz,
%      subplot(131), imagesc(EVV(:,:,k)); colorbar; impixelinfo; title('Equivolume')
%      subplot(132), imagesc(EDV(:,:,k)); colorbar; title(num2str(k)); 
%      subplot(133), imagesc(Thickness(:,:,k)); colorbar; title(num2str(k)); pause;
%  end
%  
 % Laplace schema to cover possible holes...

EVVx = ones(nx+2,ny+2,nz+2)*NaN;
EVVx(2:nx+1,2:ny+1,2:nz+1) = EVV;

EVVx(Ux==b1) = 1;
EVVx(Ux==b2) = 0;
EVVx(isnan(Ux)) = NaN;
EVVx(EVVx==0 & Ux==-Inf) = -Inf;
% 
EDVx                       = ones(nx+2,ny+2,nz+2)*NaN;
EDVx(2:nx+1,2:ny+1,2:nz+1) = EDV;

EDVx(Ux==b1)             = 1;
EDVx(Ux==b2)             = 0;
EDVx(isnan(Ux))          = NaN;
EDVx(EDVx==0 & Ux==-Inf) = -Inf;

EVV2 = laplace3d(EVVx);
EDV2 = laplace3d(EDVx);
% 
EVV2 = EVV2(2:end-1,2:end-1,2:end-1);
EDV2 = EDV2(2:end-1,2:end-1,2:end-1);