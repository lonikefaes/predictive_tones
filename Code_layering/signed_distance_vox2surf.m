function voxel_space = signed_distance_vox2surf(FV,voxel_space)

[nx, ny, nz] = size(voxel_space);


N  = patchnormals(FV);
P  = FV.vertices;

blocks   = round(linspace(1,nx*ny*nz,100));
h0       = tic;

for p = 1:length(blocks)-1,
    
    parfor k = blocks(p):blocks(p+1),
        
        voxel_space(k) = sign_dist(P,N,k,[nx,ny,nz]);
    end
    
    h1 = toc(h0);
    disp(['Finished: ',num2str(p+1),'%, ', num2str(h1/60),' minutes elapsed...']);
    
end

function signed_d = sign_dist(P,N,i,dims),

[x, y, z] = ind2sub(dims,i);
D         = sqrt((P(:,1)-(y)).^2 + (P(:,2)-(x)).^2 + (P(:,3)-(z)).^2);
[~,I]     = min(D,[],1);
p         = P(I,:);
n         = N(I,:);
% Uses Hessian form, ie : N.p = d
% I this case, d can be found as : -N.p0, when N is normalized
signed_d = -sum(bsxfun(@times, n, bsxfun(@minus,p,[y x z])), 2);