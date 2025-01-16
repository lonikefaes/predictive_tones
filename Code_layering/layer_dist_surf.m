function [FB] = layer_dist(U_out,U_outS,Mask_morph,FV)

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of streamline caluculation (Equivolume&Equidistant representation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select voxels that has to be evaluated (within GM)
UU          =  zeros(nx,ny,nz);
UU(Mask_morph==1) = 1;
UU0         = ones(nx,ny,nz);
UU0(UU0==1) = [1:length(find(UU0==1))];
UU0         = UU.*UU0;
UU(UU==1)   = [1:length(find(UU==1))];
% and get their indices and subcoordinates



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

TC = zeros(size(FV.vertices,1),1);
FB = cell(size(FV.vertices,1),1);
Fd = cell(size(FV.vertices,1),1);
thick = cell(size(FV.vertices,1),1);
thickN = cell(size(FV.vertices,1),1);
fiberN= cell(size(FV.vertices,1),1);
try
    cores = feature('numcores');
    matlabpool(cores);
end

disp('Starting parallel process calculation');
blocks = round(linspace(1,size(FV.vertices,1),100));
h0=tic;
for p = 1:length(blocks)-1,
    parfor j = blocks(p):blocks(p+1),
        xyz = FV.vertices(j,[2 1 3]);
        [FB{j}] = evaluatestreamline3Dsurf(xyz,U_out,VectorF,Roi,b1,b2,opts);
    end
    
    h1=toc(h0);
    disp(['Finished: ', num2str(100*blocks(p+1)/size(FV.vertices,1)),'%, ', num2str(h1/60),' minutes elapsed...']);
end
%  % EQUIVOLUME:
%  EVV = zeros(nx,ny,nz);
%  EVV(ind0) = EV(1:length(ind0)); 
%  % EQUIDISTANT:
%  EDV = zeros(nx,ny,nz);
%  EDV(ind0) = ED(1:length(ind0));
%  
%  Thickness       = zeros(nx,ny,nz);
%  Thickness(ind0) = TC(1:length(ind0)); 
  
