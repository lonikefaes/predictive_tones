function range = roi_dim(D,interp_factor)  



% locate region of interest:
[nx0 ny0 nz0] = size(D);
Dloc          = logical(D);
ind           = find(Dloc==1);
[ix iy iz]    = ind2sub([nx0 ny0 nz0],ind);

% lower limit has to be odd number and upper even number (assuming that anatomy was upsampled 2x)
ix_min = min(ix(:)); ix_max = max(ix(:));
iy_min = min(iy(:)); iy_max = max(iy(:));
iz_min = min(iz(:)); iz_max = max(iz(:));

if ~mod(ix_min,2), ix_min = ix_min-1; end;  if mod(ix_max,2), ix_max = ix_max-1; end; 
if ~mod(iy_min,2), iy_min = iy_min-1; end;  if mod(iy_max,2), iy_max = iy_max-1; end; 
if ~mod(iz_min,2), iz_min = iz_min-1; end;  if mod(iz_max,2), iz_max = iz_max-1; end; 

mx = [ix_min:ix_max];
my = [iy_min:iy_max];
mz = [iz_min:iz_max];

% check the size...   
nx00 = nx0/2;
ny00 = ny0/2;
nz00 = nz0/2;
mat_ind2 = reshape([1:nx00*ny00*nz00],nx00,ny00,nz00);
MAT_ind2 = [];
for i = 1:nz00
    MAT_ind2 = cat(3,MAT_ind2,repmat(kron(mat_ind2(:,:,i),ones(2)),[1,1,2]));
end

Original = MAT_ind2(mx,my,mz);
indorg   = unique(Original(:));
Org      = zeros(nx00,ny00,nz00);
Org(indorg) = 1;
[xi, yi, zi] = ind2sub([nx00 ny00 nz00],indorg);

mx0 = min(xi):max(xi); % 80:115
my0 = min(yi):max(yi); % 123:165
mz0 = min(zi):max(zi); % 31:43


nx = nx0*(interp_factor/2);
ny = ny0*(interp_factor/2);
nz = nz0*(interp_factor/2); 


Up2 = [];
for i = 1:nz00
    Up2 = cat(3,Up2,repmat(kron(Org(:,:,i),ones(2)),[1,1,2]));
end

Up5 = [];
for i = 1:nz00
    Up5 = cat(3,Up5,repmat(kron(Org(:,:,i),ones(interp_factor)),[1,1,interp_factor]));
end

ind2 = find(Up2==1);

[xi2, yi2, zi2] = ind2sub([nx0 ny0 nz0],ind2);

mx2 = min(xi2):max(xi2); % 159:230;
my2 = min(yi2):max(yi2); % 245:330;
mz2 = min(zi2):max(zi2); % 61:86

ind5 = find(Up5==1);

[xi5, yi5, zi5] = ind2sub([nx ny nz],ind5);

mx5 = min(xi5):max(xi5); % 317:460
my5 = min(yi5):max(yi5); % 489:660
mz5 = min(zi5):max(zi5); % 121:172

range.r0.mx = mx0; 
range.r0.my = my0; 
range.r0.mz = mz0; 
range.r2.mx = mx2; 
range.r2.my = my2; 
range.r2.mz = mz2; 
range.r5.mx = mx5; 
range.r5.my = my5; 
range.r5.mz = mz5;

save([cd,'\FCN\range.mat'],'range');