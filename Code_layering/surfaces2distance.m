function [DM TT CC] = surfaces2distance(SP,DV,Mask,Thickness,Curvature,SP2),

[nx, ny, nz] = size(Mask);

Points0  = cell2mat(SP');
Points1  = cell2mat(SP2');
D       = cell2mat(DV');
T       = cell2mat(Thickness)';
C       = cell2mat(Curvature)';
WM      = Mask==2;
GM0     = Mask==1;
CSF     = Mask==0;

[xx,yy,zz] = ndgrid(-2:2);
nhood      = sqrt(xx.^2 + yy.^2 + zz.^2) <= 2;

GM = imdilate(GM0,nhood);

WM  = WM.*~GM;
CSF = CSF.*~GM;

IndWM               = find(WM==1);
IndCSF              = find(CSF==1);

[wm_y,wm_x,wm_z]    = ind2sub([nx,ny,nz],IndWM);
[csf_y,csf_x,csf_z] = ind2sub([nx,ny,nz],IndCSF);

PointsWM            = [wm_x,wm_y,wm_z];
PointsCSF           = [csf_x,csf_y,csf_z];

Dwm                 = -0.01*ones(size(PointsWM,1),1);
Dcsf                =  1.01*ones(size(PointsCSF,1),1);


[X,Y,Z] = meshgrid([1:ny],[1:nx],[1:nz]);

Points  = [Points0;PointsWM;PointsCSF];
D       = [D;Dwm;Dcsf];

DM      = griddata(Points(:,1),Points(:,2),Points(:,3),D,X,Y,Z,'linear');

DM(isnan(DM))  = 0;
DM([1,nx],:,:) = NaN;
DM(:,[1,ny],:) = NaN;
DM(:,:,[1,nz]) = NaN;


TT = zeros(nx,ny,nz);
TTn = zeros(nx,ny,nz);

CC = zeros(nx,ny,nz);
CCn = zeros(nx,ny,nz);
for i = 1:length(T)
    
   TT(round(Points0(i,2)),round(Points0(i,1)),round(Points0(i,3))) = TT(round(Points0(i,2)),round(Points0(i,1)),round(Points0(i,3))) + T(i); 
   TTn(round(Points0(i,2)),round(Points0(i,1)),round(Points0(i,3))) = TTn(round(Points0(i,2)),round(Points0(i,1)),round(Points0(i,3))) + 1; 

end

for i = 1:length(C)
  CC(round(Points1(i,2)),round(Points1(i,1)),round(Points1(i,3))) = CC(round(Points1(i,2)),round(Points1(i,1)),round(Points1(i,3))) + C(i); 
  CCn(round(Points1(i,2)),round(Points1(i,1)),round(Points1(i,3))) = CCn(round(Points1(i,2)),round(Points1(i,1)),round(Points1(i,3))) + 1;
end
TT = TT./TTn;
TT(isnan(TT)) = 0;

CC            = CC./CCn;
CC(isnan(CC)) = 0;
%laplace3d()

